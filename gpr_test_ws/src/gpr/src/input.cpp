#include "input.h"

void Imu_Data_t::feed(sensor_msgs::ImuConstPtr pMsg)
{
    ros::Time now = ros::Time::now();

    msg = *pMsg;
    rcv_stamp = now;

    w(0) = msg.angular_velocity.x;
    w(1) = msg.angular_velocity.y;
    w(2) = msg.angular_velocity.z;

    a(0) = msg.linear_acceleration.x;
    a(1) = msg.linear_acceleration.y;
    a(2) = msg.linear_acceleration.z;

    q.x() = msg.orientation.x;
    q.y() = msg.orientation.y;
    q.z() = msg.orientation.z;
    q.w() = msg.orientation.w;

    // check the frequency
    static int one_min_count = 9999;
    static ros::Time last_clear_count_time = ros::Time(0.0);
    if ( (now - last_clear_count_time).toSec() > 1.0 )
    {
        if ( one_min_count < 100 )
        {
            ROS_WARN("IMU frequency seems lower than 100Hz, which is too low!");
        }
        one_min_count = 0;
        last_clear_count_time = now;
    }
    one_min_count ++;
}

Odom_Data_t::Odom_Data_t()
{
    rcv_stamp = ros::Time(0);
    q.setIdentity();
    recv_new_msg = false;
};

void Odom_Data_t::feed(nav_msgs::OdometryConstPtr pMsg)
{
    ros::Time now = ros::Time::now();

    msg = *pMsg;
    rcv_stamp = now;
    recv_new_msg = true;

    p(0) = msg.pose.pose.position.x;
    p(1) = msg.pose.pose.position.y;
    p(2) = msg.pose.pose.position.z;

    v(0) = msg.twist.twist.linear.x;
    v(1) = msg.twist.twist.linear.y;
    v(2) = msg.twist.twist.linear.z;

    q.x() = msg.pose.pose.orientation.x;
    q.y() = msg.pose.pose.orientation.y;
    q.z() = msg.pose.pose.orientation.z;
    q.w() = msg.pose.pose.orientation.w;

// #define VEL_IN_BODY
#ifdef VEL_IN_BODY /* Set to 1 if the velocity in odom topic is relative to current body frame, not to world frame.*/
    Eigen::Quaternion<double> wRb_q(msg.pose.pose.orientation.w, msg.pose.pose.orientation.x, msg.pose.pose.orientation.y, msg.pose.pose.orientation.z);
    Eigen::Matrix3d wRb = wRb_q.matrix();
    v = wRb * v;

    static int count = 0;
    if (count++ % 500 == 0)
        ROS_WARN("VEL_IN_BODY!!!");
#endif

    // check the frequency
    static int one_min_count = 9999;
    static ros::Time last_clear_count_time = ros::Time(0.0);
    if ( (now - last_clear_count_time).toSec() > 1.0 )
    {
        if ( one_min_count < 100 )
        {
            ROS_WARN("ODOM frequency seems lower than 100Hz, which is too low!");
        }
        one_min_count = 0;
        last_clear_count_time = now;
    }
    one_min_count ++;
}

void RC_Data_t::feed(mavros_msgs::RCInConstPtr pMsg)
{
    msg = *pMsg;
    rcv_stamp = ros::Time::now();

    double channel_7 = ((double)msg.channels[6] - 1500.0);
    if (channel_7 < 0)
    {
        est_trigger = true;
    }
    else if (channel_7 > 0 || channel_7 == 0)
    {
        est_trigger = false;
    }
}

void GPS_VEL_Data_t::feed(geometry_msgs::TwistStampedConstPtr pMsg)
{
    ros::Time now = ros::Time::now();
    
    msg = *pMsg;
    rcv_stamp = now;

    gps_vel(0) = msg.twist.linear.x;
    gps_vel(1) = msg.twist.linear.y;
    gps_vel(2) = msg.twist.linear.z;
}

void Mag_Data_t::feed(sensor_msgs::MagneticFieldConstPtr pMsg)
{
    ros::Time now = ros::Time::now();
    
    msg = *pMsg;
    rcv_stamp = now;

    mag(0) = msg.magnetic_field.x;
    mag(1) = msg.magnetic_field.y;
    mag(2) = msg.magnetic_field.z;
}