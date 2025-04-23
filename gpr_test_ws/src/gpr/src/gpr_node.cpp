#include <iostream>
#include <ros/ros.h>

#include "gpr.h"

using namespace std;

vector<double> manual_correlate(const vector<double>& a, const vector<double>& v) {
    vector<double> a_copy = a;
    vector<double> v_copy = v;

    // Ensure v is the shorter sequence
    if (v_copy.size() > a_copy.size()) {
        swap(a_copy, v_copy);
    }

    int a_len = a_copy.size();
    int v_len = v_copy.size();
    int valid_len = a_len + v_len - 1;
    
    vector<double> result(valid_len, 0.0);

    // Flip v for correlation
    vector<double> v_flipped(v_len);
    for (int i = 0; i < v_len; i++) {
        v_flipped[i] = v_copy[v_len - 1 - i];
    }

    // Compute valid correlation
    for (int i = 0; i < valid_len; i++) {
        for (int j = 0; j < v_len; j++) {
            if (0 <= i - j && i - j < a_len){
                result[i] += a_copy[i - j] * v_flipped[j];
            }
        }
    }

    return result;
}

int main(int argc,char **argv)
{
    ros::init(argc,argv,"gpr");
    ros::NodeHandle n("~");
    ros::Time last_rcv_stamp = ros::Time::now();
    double init_rcv_stamp = last_rcv_stamp.toSec();

    Imu_Data_t imu_data;
    RC_Data_t rc_data;
    GPS_VEL_Data_t gps_vel_data;
    Mag_Data_t mag_data;

    Butterworth_lowpass_filter acc_x_blf, acc_y_blf, acc_z_blf, agv_x_blf, agv_y_blf, agv_z_blf;
    Butterworth_lowpass_filter gps_vel_x_blf, gps_vel_y_blf, gps_vel_z_blf;
    Freq_estimator freq_est;
    Phase_based_GPR p_GPR;
    EKF_marg ekf_estimator;
    ekf_estimator.EKF_init();

    // mavros topic defined
    ros::Subscriber imu_sub = n.subscribe<sensor_msgs::Imu>("/mavros/imu/data_raw", 100, 
                                                             boost::bind(&Imu_Data_t::feed, &imu_data, _1),
                                                             ros::VoidConstPtr(),
                                                             ros::TransportHints().tcpNoDelay());
    ros::Subscriber rc_sub = n.subscribe<mavros_msgs::RCIn>("/mavros/rc/in",
                                                              10,
                                                              boost::bind(&RC_Data_t::feed, &rc_data, _1));
    ros::Subscriber gps_vel_sub = n.subscribe<geometry_msgs::TwistStamped>("/mavros/global_position/raw/gps_vel", 100,
                                                                          boost::bind(&GPS_VEL_Data_t::feed, &gps_vel_data, _1),
                                                                          ros::VoidConstPtr(),
                                                                          ros::TransportHints().tcpNoDelay());
    ros::Subscriber mag_sub = n.subscribe<sensor_msgs::MagneticField>("/mavros/imu/mag", 100,
                                                                     boost::bind(&Mag_Data_t::feed, &mag_data, _1),
                                                                     ros::VoidConstPtr(),
                                                                     ros::TransportHints().tcpNoDelay());
    ros::Publisher imu_osc_pub = n.advertise<sensor_msgs::Imu>("/imu_osc", 100);
    ros::Publisher imu_no_osc_pub = n.advertise<sensor_msgs::Imu>("/imu_no_osc", 100);

    vector<double> time_list, acc_z_list, agv_y_list;
    vector<double> time_list_test_corr, signal_est, signal_raw, lag_max_list;
    double phase_curr = 0.0, lag_max = 0.0;
    vector<Phase_based_data> p_data_list;

    double time_gps_last = 0, gps_vel_x_last = 0, gps_vel_y_last = 0, gps_vel_z_last = 0, time_last_ekf = 0;
    Eigen::Vector3d gps_acc;
    double gps_acc_x = 0, gps_acc_y = 0, gps_acc_z = 0;
    
    ros::Rate rate(200.0);
    while(ros::ok())
    {
        rate.sleep();
        ros::spinOnce();
        
        // imu pub
        sensor_msgs::Imu Imu_osc;
        sensor_msgs::Imu Imu_no_osc;

        time_list.push_back(imu_data.rcv_stamp.toSec() - init_rcv_stamp);
        acc_z_list.push_back(imu_data.a(2));
        agv_y_list.push_back(imu_data.w(1));
        if (time_list.size() > 512)
        {
            time_list.erase(time_list.begin());
            acc_z_list.erase(acc_z_list.begin());
            agv_y_list.erase(agv_y_list.begin());
        }

        if (gps_vel_data.rcv_stamp.toSec() != time_gps_last)
        {
            double dt_gps = gps_vel_data.rcv_stamp.toSec() - time_gps_last;
            double gps_vel_x_filtered = gps_vel_x_blf.process_sample(gps_vel_data.gps_vel(0), 10, 100);
            double gps_vel_y_filtered = gps_vel_y_blf.process_sample(gps_vel_data.gps_vel(1), 10, 100);
            double gps_vel_z_filtered = gps_vel_z_blf.process_sample(gps_vel_data.gps_vel(2), 10, 100);
            gps_acc_x = (gps_vel_x_filtered - gps_vel_x_last) / (dt_gps + 1e-8);
            gps_acc_y = (gps_vel_y_filtered - gps_vel_y_last) / (dt_gps + 1e-8);
            gps_acc_z = (gps_vel_z_filtered - gps_vel_z_last) / (dt_gps + 1e-8);
            gps_vel_x_last = gps_vel_x_filtered;
            gps_vel_y_last = gps_vel_y_filtered;
            gps_vel_z_last = gps_vel_z_filtered;
            time_gps_last = gps_vel_data.rcv_stamp.toSec();
        }
        gps_acc << gps_acc_x, gps_acc_y, gps_acc_z;

        if (!rc_data.est_trigger)
        {
            ROS_WARN("START ESTIMATION!");
            // calculate frequency
            struct Freq_est_data freq_data;
            freq_data = freq_est.freq_estimator(time_list, acc_z_list, agv_y_list);

            // calculate current phase
            phase_curr += (time_list[time_list.size() - 1] - time_list[time_list.size() - 2]) * 2 * M_PI * freq_data.freq_mean_4_freq;
            phase_curr += max(min((time_list[511] - time_list[510]) * 2 * M_PI * freq_data.freq_mean_4_freq * 0.1 * lag_max, M_PI/10), -M_PI/10);
            phase_curr = fmod(phase_curr, 2 * M_PI);

            struct Phase_based_data p_data;
            p_data.time_record = imu_data.rcv_stamp.toSec() - init_rcv_stamp;
            p_data.phase = phase_curr;
            p_data.acc_x_imu_data = imu_data.a(0); 
            p_data.acc_y_imu_data = imu_data.a(1); 
            p_data.acc_z_imu_data = imu_data.a(2); 
            p_data.agv_x_imu_data = imu_data.w(0);
            p_data.agv_y_imu_data = imu_data.w(1);
            p_data.agv_z_imu_data = imu_data.w(2);
            p_data_list.push_back(p_data);

            if (p_data_list.size() > 224)
            {
                p_data_list.erase(p_data_list.begin());
                
                // Gassuain Process Regression
                struct f_est_var osc_est_var;
                osc_est_var = p_GPR.update_estimation_essentials(p_data_list);

                // publish imu osc message
                Imu_osc.header.stamp = ros::Time::now();
                Imu_osc.header.frame_id = "imu_frame";
                Imu_osc.linear_acceleration.x = osc_est_var.acc_osc(0);
                Imu_osc.linear_acceleration.y = osc_est_var.acc_osc(1);
                Imu_osc.linear_acceleration.z = osc_est_var.acc_osc(2);
                Imu_osc.angular_velocity.x = osc_est_var.agv_osc(0);
                Imu_osc.angular_velocity.y = osc_est_var.agv_osc(1);
                Imu_osc.angular_velocity.z = osc_est_var.agv_osc(2);
                imu_osc_pub.publish(Imu_osc);

                // butterworth filter for no osc data
                double acc_x_no_osc = imu_data.a(0) - osc_est_var.acc_osc(0);
                double acc_x_no_osc_filtered = acc_x_blf.process_sample(acc_x_no_osc, 10, 200);
                double acc_y_no_osc = imu_data.a(1) - osc_est_var.acc_osc(1);
                double acc_y_no_osc_filtered = acc_y_blf.process_sample(acc_y_no_osc, 10, 200);
                double acc_z_no_osc = imu_data.a(2) - osc_est_var.acc_osc(2);
                double acc_z_no_osc_filtered = acc_z_blf.process_sample(acc_z_no_osc, 10, 200);
                double agv_x_no_osc = imu_data.w(0) - osc_est_var.agv_osc(0);
                double agv_x_no_osc_filtered = agv_x_blf.process_sample(agv_x_no_osc, 10, 200);
                double agv_y_no_osc = imu_data.w(1) - osc_est_var.agv_osc(1);
                double agv_y_no_osc_filtered = agv_y_blf.process_sample(agv_y_no_osc, 10, 200);
                double agv_z_no_osc = imu_data.w(2) - osc_est_var.agv_osc(2);
                double agv_z_no_osc_filtered = agv_z_blf.process_sample(agv_z_no_osc, 10, 200);
                Eigen::Vector3d acc_no_osc, agv_no_osc;
                acc_no_osc << acc_x_no_osc_filtered, acc_y_no_osc_filtered, acc_z_no_osc_filtered;
                agv_no_osc << agv_x_no_osc_filtered, agv_y_no_osc_filtered, agv_z_no_osc_filtered;

                // calculate correlation
                signal_est.push_back(osc_est_var.acc_osc(2));
                signal_raw.push_back(p_data_list[p_data_list.size() - 1].acc_z_imu_data);
                if (signal_est.size() > 40)
                {
                    signal_est.erase(signal_est.begin());
                    signal_raw.erase(signal_raw.begin());
                    vector<double> corr = manual_correlate(signal_est, signal_raw);
                    vector<double> lags;
                    for (int i = -39; i < 40; i++)
                    {
                        lags.push_back(i);
                    }
                    int corr_size = corr.size();
                    auto maxIt = max_element(corr.begin(), corr.end());
                    int maxcorr_Index = distance(corr.begin(), maxIt);
                    double lag_max_candidate = lags[maxcorr_Index];
                    if (lag_max_list.size() > 0 && abs(lag_max_candidate - lag_max_list[lag_max_list.size() - 1]) > 10)
                    {
                        lag_max = lag_max_list[lag_max_list.size() - 1];
                    }
                    else
                    {
                        lag_max = lag_max_candidate;
                    }
                    lag_max_list.push_back(lag_max);
                }

                // quaternion estimate by EKF
                double mag_x_norm = mag_data.mag(0) / sqrt(mag_data.mag(0) * mag_data.mag(0) + mag_data.mag(1) * mag_data.mag(1) + mag_data.mag(2) * mag_data.mag(2));
                double mag_y_norm = mag_data.mag(1) / sqrt(mag_data.mag(0) * mag_data.mag(0) + mag_data.mag(1) * mag_data.mag(1) + mag_data.mag(2) * mag_data.mag(2));
                double mag_z_norm = mag_data.mag(2) / sqrt(mag_data.mag(0) * mag_data.mag(0) + mag_data.mag(1) * mag_data.mag(1) + mag_data.mag(2) * mag_data.mag(2));
                Eigen::Vector3d mag_norm;
                mag_norm << mag_x_norm, mag_y_norm, mag_z_norm;
                Eigen::Quaterniond quat_estimated;
                quat_estimated.w() = ekf_estimator.x(0);
                quat_estimated.x() = ekf_estimator.x(1);
                quat_estimated.y() = ekf_estimator.x(2);
                quat_estimated.z() = ekf_estimator.x(3);
                Eigen::Matrix3d R_estimated;
                R_estimated = quat_estimated.toRotationMatrix();
                Eigen::Vector3d acc_gps = R_estimated.transpose() * gps_acc;
                double acc_x_rect = acc_x_no_osc_filtered - acc_gps(0);
                double acc_y_rect = acc_y_no_osc_filtered - acc_gps(1);
                double acc_z_rect = acc_z_no_osc_filtered - acc_gps(2);
                double acc_x_rect_norm = acc_x_rect / sqrt(acc_x_rect * acc_x_rect + acc_y_rect * acc_y_rect + acc_z_rect * acc_z_rect);
                double acc_y_rect_norm = acc_y_rect / sqrt(acc_x_rect * acc_x_rect + acc_y_rect * acc_y_rect + acc_z_rect * acc_z_rect);
                double acc_z_rect_norm = acc_z_rect / sqrt(acc_x_rect * acc_x_rect + acc_y_rect * acc_y_rect + acc_z_rect * acc_z_rect);
                Eigen::Vector3d acc_rect_norm;
                acc_rect_norm << acc_x_rect_norm, acc_y_rect_norm, acc_z_rect_norm;
                ros::Time ekf_time = ros::Time::now();
                double dt_ekf = ekf_time.toSec() - time_last_ekf;
                time_last_ekf = ekf_time.toSec();
                Eigen::Vector3d mag_new;
                mag_new = R_estimated * mag_norm;
                mag_new(2) = 0.0;
                double mag_new_norm = sqrt(mag_new(0) * mag_new(0) + mag_new(1) * mag_new(1) + mag_new(2) * mag_new(2));
                Eigen::Vector3d mag_new_proj;
                mag_new_proj << mag_new(0) / mag_new_norm, mag_new(1) / mag_new_norm, mag_new(2) / mag_new_norm;
                Eigen::Vector3d mag_new_proj_back;
                mag_new_proj_back = R_estimated.transpose() * mag_new_proj;
                ekf_estimator.predict(agv_no_osc, dt_ekf);
                ekf_estimator.update(acc_rect_norm, mag_new_proj_back);

                // publish cycle average imu message
                Imu_no_osc.header.stamp = ros::Time::now();
                Imu_no_osc.header.frame_id = "imu_frame";
                Imu_no_osc.linear_acceleration.x = acc_x_no_osc_filtered;
                Imu_no_osc.linear_acceleration.y = acc_y_no_osc_filtered;
                Imu_no_osc.linear_acceleration.z = acc_z_no_osc_filtered;
                Imu_no_osc.orientation.w = ekf_estimator.x(0);
                Imu_no_osc.orientation.x = ekf_estimator.x(1);
                Imu_no_osc.orientation.y = ekf_estimator.x(2);
                Imu_no_osc.orientation.z = ekf_estimator.x(3);
                Imu_no_osc.angular_velocity.x = agv_x_no_osc_filtered;
                Imu_no_osc.angular_velocity.y = agv_y_no_osc_filtered;
                Imu_no_osc.angular_velocity.z = agv_z_no_osc_filtered;
                imu_no_osc_pub.publish(Imu_no_osc);
            }
        }
    }

    return 0;
}