#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <complex>
#include <random>

#include "input.h"

using namespace std;

class Butterworth_lowpass_filter
{
public:
    double x[3] = {0, 0, 0};
    double y[2] = {0, 0};

    double process_sample(double sample, double cutoff_freq, double sample_freq);
};

struct Freq_est_data
{
    double time_4_freq;
    double freq_mean_4_freq;
    double freq_ob_noise_4_freq;
};

class Freq_estimator
{
public:
    double start_freq = 1.0;
    double end_freq = 8.0;
    double ignore_mag_norm = 0.2;

    vector<double> fftfreq_manual(const double &n, const double &d);
    vector<complex<double>> radix2_fft(const vector<complex<double>>& x);
    Freq_est_data freq_estimator(const vector<double> &time_list, const vector<double> &acc_z_list, const vector<double> &arg_y_list);
};

struct mean_var
{
    double mean;
    double var;
};

struct k_means_return
{
    vector<double> cluster_centers;
    vector<int> labels;
};

struct Phase_based_data
{
    double time_record;
    double phase;
    double acc_x_imu_data;
    double acc_y_imu_data;
    double acc_z_imu_data;
    double agv_x_imu_data;
    double agv_y_imu_data;
    double agv_z_imu_data;
    double sigma_acc_x;
    double sigma_acc_y;
    double sigma_acc_z;
    double sigma_agv_x;
    double sigma_agv_y;
    double sigma_agv_z;
};

struct f_est_var
{
    Eigen::Vector3d acc_osc;
    Eigen::Vector3d agv_osc;
    Eigen::Vector3d acc_var;
    Eigen::Vector3d agv_var;
};

class Phase_based_GPR
{
public:
    // params for k_means
    int cluster_num = 15.0;
    double Fourier_kernel_trunc_order = 6.0;

    vector<double> phase_output;
    vector<double> cluster_centers;
    bool cluster = false;

    double vector_mean(const vector<double> &a);
    mean_var mean_variance(const vector<double> &a);
    vector<int> sort_indexes(vector< double> v);
    double similarity_defined_by_chainplain(const double &phase, const double &cluster_centers);
    k_means_return phase_kmeans_fast(const vector<double> &phase);
    double Fourier_kernel(const double &phi_a, const double &phi_b);
    f_est_var update_estimation_essentials(const vector<Phase_based_data> &imu_data_input);
};

class EKF_marg
{
public:
    Eigen::VectorXd x;
    Eigen::MatrixXd P;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd R;

    Eigen::Vector3d gravitational_field;
    Eigen::Vector3d magnetic_field;

    void EKF_init();
    Eigen::VectorXd quaternion_multiply(const Eigen::VectorXd &q1, const Eigen::VectorXd &q2);
    void predict(const Eigen::Vector3d &gyro_data, const double &dt);
    void update(const Eigen::Vector3d &acc_data, const Eigen::Vector3d &mag_data);
};