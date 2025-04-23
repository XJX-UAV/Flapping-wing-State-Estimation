#include "gpr.h"

double Butterworth_lowpass_filter::process_sample(double sample, double cutoff_freq, double sample_freq)
{
    double fr = sample_freq / cutoff_freq;
    double ohm = tan(M_PI / fr);
    double c = 1.0 + 1.414 * ohm + ohm * ohm;

    double b_a[2][3];
    b_a[0][0] = ohm * ohm / c;
    b_a[0][1] = 2 * b_a[0][0];
    b_a[0][2] = b_a[0][0];
    b_a[1][0] = 1;
    b_a[1][1] = 2 * (ohm*ohm - 1) / c;
    b_a[1][2] = (1 - 1.414 * ohm + ohm * ohm) / c;
    
    x[2] = x[1];
    x[1] = x[0];
    x[0] = sample;

    double y_new = ((b_a[0][0] * x[0] + b_a[0][1] * x[1] + b_a[0][2] * x[2]) - (b_a[1][1] * y[0] + b_a[1][2] * y[1]));
    y[1] = y[0];
    y[0] = y_new;

    return y_new;
}

vector<double> Freq_estimator::fftfreq_manual(const double &N, const double &d)
{
    vector<double> result;
    for (int i = 0; i < N / 2; i++)
    {
        result.push_back(i / (N * d));
    }

    return result;
}

vector<complex<double>> Freq_estimator::radix2_fft(const vector<complex<double>>& x) {
    int N = x.size();
    if (N <= 1) return x;
    
    vector<complex<double>> even, odd;
    for (int i = 0; i < N; i += 2) {
        even.push_back(x[i]);
        if (i + 1 < N) odd.push_back(x[i + 1]);
    }
    
    vector<complex<double>> even_fft = radix2_fft(even);
    vector<complex<double>> odd_fft = radix2_fft(odd);
    
    vector<complex<double>> result(N);
    for (int k = 0; k < N / 2; ++k) {
        complex<double> twiddle = exp(complex<double>(0, -2.0 * M_PI * k / N)) * odd_fft[k];
        result[k] = even_fft[k] + twiddle;
        result[k + N / 2] = even_fft[k] - twiddle;
    }
    return result;
}

Freq_est_data Freq_estimator::freq_estimator(const vector<double> &time_list, const vector<double> &acc_z_list, const vector<double> &agv_y_list)
{
    struct Freq_est_data freq_data;
    // trasform data type from double to complex
    vector<complex<double>> acc_z_list_complex, agv_y_list_complex;
    for (int i = 0; i < acc_z_list.size(); i++)
    {
        complex<double> acc_z_complex = acc_z_list[i];
        acc_z_list_complex.push_back(acc_z_complex);
    }
    for (int i = 0; i < agv_y_list.size(); i++)
    {
        complex<double> arg_y_complex = agv_y_list[i];
        agv_y_list_complex.push_back(arg_y_complex);
    }

    // fft
    double inter = (time_list[time_list.size() - 1] - time_list[0]) / 512;
    vector<double> acc_z_fft_result, acc_z_fft_result_new, agv_y_fft_result, agv_y_fft_result_new;
    vector<int> iter_list;
    vector<double> freqs_raw, freqs;
    freqs_raw = fftfreq_manual(time_list.size(), inter);
    int iter = 0;
    while (freqs_raw[iter] < end_freq)
    {
        if (freqs_raw[iter] > start_freq)
        {
            freqs.push_back(freqs_raw[iter]);
            iter_list.push_back(iter);
        }
        iter += 1;
    }
    int i_start = iter_list[0];
    int i_end = iter_list.back();
    // fft acc
    vector<complex<double>> acc_z_fft_result_raw = radix2_fft(acc_z_list_complex);
    for (int j = 0; j < acc_z_fft_result_raw.size() / 2; j++)
    {
        acc_z_fft_result.push_back(norm(acc_z_fft_result_raw[j]));
    }

    for (int i = i_start; i < i_end + 1; i++)
    {
        double acc_z_fft_result_inter = acc_z_fft_result[i];
        acc_z_fft_result_new.push_back(acc_z_fft_result_inter);
    }
    double acc_z_max = *max_element(acc_z_fft_result_new.begin(), acc_z_fft_result_new.end());
    for(int j = 0; j < acc_z_fft_result_new.size(); j++)
    {
        acc_z_fft_result_new[j] = acc_z_fft_result_new[j] / acc_z_max;
        if (acc_z_fft_result_new[j] < ignore_mag_norm)
        {
            acc_z_fft_result_new[j] = 0.01;
        }
    }
    double acc_z_sum = 0;
    for (int k = 0; k < acc_z_fft_result_new.size(); k++)
    {
        acc_z_sum += acc_z_fft_result_new[k];
    }
    double freq_mean_acc_z_sum = 0;
    for (int m = 0; m < freqs.size(); m++)
    {
        freq_mean_acc_z_sum += freqs[m] * acc_z_fft_result_new[m];
    }
    double freq_mean_acc_z = freq_mean_acc_z_sum / acc_z_sum;

    vector<double> freq_error_acc_z;
    for (int i = 0; i < freqs.size(); i++)
    {
        double freq_error = (freqs[i] - freq_mean_acc_z) * (freqs[i] - freq_mean_acc_z);
        freq_error_acc_z.push_back(freq_error);
    }
    double freq_var_acc_z_sum = 0;
    for (int n = 0; n < freq_error_acc_z.size(); n++)
    {
        freq_var_acc_z_sum += freq_error_acc_z[n] * acc_z_fft_result_new[n];
    }
    double freq_var_acc_z = freq_var_acc_z_sum / freq_mean_acc_z_sum;

    // fft agv
    vector<complex<double>> agv_y_fft_result_raw = radix2_fft(agv_y_list_complex);
    for (int j = 0; j < agv_y_fft_result_raw.size() / 2; j++)
    {
        agv_y_fft_result.push_back(norm(agv_y_fft_result_raw[j]));
    }

    for (int i = i_start; i < i_end + 1; i++)
    {
        double agv_y_fft_result_inter = agv_y_fft_result[i];
        agv_y_fft_result_new.push_back(agv_y_fft_result_inter);
    }

    double agv_y_max = *max_element(agv_y_fft_result_new.begin(), agv_y_fft_result_new.end());
    for(int j = 0; j < agv_y_fft_result_new.size(); j++)
    {
        agv_y_fft_result_new[j] = agv_y_fft_result_new[j] / agv_y_max;
        if (agv_y_fft_result_new[j] < ignore_mag_norm)
        {
            agv_y_fft_result_new[j] = 0.01;
        }
    }
    double agv_y_sum = 0;
    for (int k = 0; k < agv_y_fft_result_new.size(); k++)
    {
        agv_y_sum += agv_y_fft_result_new[k];
    }
    double freq_mean_agv_y_sum = 0;
    for (int m = 0; m < freqs.size(); m++)
    {
        freq_mean_agv_y_sum += freqs[m] * agv_y_fft_result_new[m];
    }
    double freq_mean_agv_y = freq_mean_agv_y_sum / agv_y_sum;
    vector<double> freq_error_agv_y;
    for (int i = 0; i < freqs.size(); i++)
    {
        double freq_error = (freqs[i] - freq_mean_agv_y) * (freqs[i] - freq_mean_agv_y);
        freq_error_agv_y.push_back(freq_error);
    }
    double freq_var_agv_y_sum = 0;
    for (int n = 0; n < freq_error_agv_y.size(); n++)
    {
        freq_var_agv_y_sum += freq_error_agv_y[n] * agv_y_fft_result_new[n];
    }
    double freq_var_agv_y = freq_var_agv_y_sum / freq_mean_agv_y_sum;

    // fuse acc_z and agv_y
    double freq_mean_est = (freq_mean_acc_z * freq_var_acc_z + freq_mean_agv_y * freq_var_agv_y) / (freq_var_acc_z + freq_var_agv_y);
    double freq_var_est = (freq_var_acc_z * freq_var_acc_z + freq_var_agv_y * freq_var_agv_y) / (freq_var_acc_z + freq_var_agv_y);
    double len = time_list.size();

    freq_data.time_4_freq = time_list[len - 1];
    freq_data.freq_mean_4_freq = freq_mean_est;
    freq_data.freq_ob_noise_4_freq = freq_var_est;
    return freq_data;
}

double Phase_based_GPR::vector_mean(const vector<double> &a)
{
    double sum_ = 0.0;
    for (int i = 0; i < a.size(); i++)
    {
        sum_ += a[i];
    }
    sum_ = sum_ / a.size();
    return sum_;
}

mean_var Phase_based_GPR::mean_variance(const vector<double> &a)
{
    struct mean_var A;
    double mean = 0.0;
    for (int j = 0; j < a.size(); j++)
    {
        mean += a[j];
    }
    mean = mean / a.size();
    double variance = 0.0;
    for (int i = 0 ; i < a.size() ; i++)
    {
        variance = variance + (a[i]- mean) * (a[i] - mean);
    }
    variance = variance / a.size();
    double standard_deviation = sqrt(variance);
    A.mean = mean;
    A.var = standard_deviation;
    return A;
}

vector<int> Phase_based_GPR::sort_indexes(vector<double> v) {

  // initialize original index locations
  vector<int> idx(v.size());
  for (int i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [& v](size_t i1, size_t i2) {return v[i1] <  v[i2];});

  return idx;
}

double Phase_based_GPR::similarity_defined_by_chainplain(const double &phase, const double &cluster_centers)
{
    double similarity = 1 - cos((phase - cluster_centers) / 2.0);
    return similarity;
}

k_means_return Phase_based_GPR::phase_kmeans_fast(const vector<double> &phase)
{
    struct k_means_return k_means;
    // Randomly initialize cluster centers (use phase values as initial centers)
    vector<int> max_simi_labels;
    vector<double> phase_random = phase;
    if (cluster == false)
    {
        srand(unsigned(time(NULL)));
        random_shuffle(phase_random.begin(), phase_random.end());
        for (int i = 0; i < cluster_num; i++)
        {
            cluster_centers.push_back(phase_random[i]);
        }
        cluster = true;
    }

    int iteration = 0;
    while (iteration < 100)
    {
        double similarity[phase.size()][cluster_num];
        for (int i = 0; i < phase.size(); i++)
        {
            for (int j = 0; j < int(cluster_num); j++)
            {
                similarity[i][j] = similarity_defined_by_chainplain(phase[i], cluster_centers[j]);
            }
        }
        
        for (int k = 0; k < phase.size(); k++)
        {
            int minPosition = min_element(similarity[k] ,similarity[k] + cluster_num) - similarity[k];
            max_simi_labels.push_back(minPosition);
        }
        vector<double> new_cluster_centers;
        for (int l = 0; l < cluster_num; l++)
        {
            vector<double> phase_label;
            for (int m = 0; m < phase.size(); m++)
            {
                if (max_simi_labels[m] == l)
                {
                    phase_label.push_back(phase[m]);
                }
            }
            double mean = vector_mean(phase_label);
            new_cluster_centers.push_back(mean);
        }

        bool convergence = true;
        for (int n = 0; n < cluster_num; n++)
        {
            if (new_cluster_centers[n] - cluster_centers[n] > 0.0001)
            {
                convergence = false;
            }
        }
        if (convergence == true)
        {
            break;
        }

        cluster_centers = new_cluster_centers;

        iteration += 1;
    }
    vector<int> max_simi_labels_new;
    vector<int> sorted_indices;
    sorted_indices = sort_indexes(cluster_centers);
    sort(cluster_centers.begin(), cluster_centers.end());
    for (int i = 0; i < max_simi_labels.size(); i++)
    {
        max_simi_labels_new.push_back(sorted_indices[max_simi_labels[i]]);
    }
    
    k_means.cluster_centers = cluster_centers;
    k_means.labels = max_simi_labels_new;
    return k_means;
}

double Phase_based_GPR::Fourier_kernel(const double &phi_a, const double &phi_b)
{
    double fourier_vec_a[int(2 * Fourier_kernel_trunc_order)];
    double fourier_vec_b[int(2 * Fourier_kernel_trunc_order)];
    for (int i = 1; i < Fourier_kernel_trunc_order + 1; i++)
    {
        fourier_vec_a[2 * i - 2] = sin(i * phi_a);
        fourier_vec_a[2 * i - 1] = cos(i * phi_a);
        fourier_vec_b[2 * i - 2] = sin(i * phi_b);
        fourier_vec_b[2 * i - 1] = cos(i * phi_b);
    }
    
    double raw_output = 0.0;
    for (int j = 0; j < 2 * Fourier_kernel_trunc_order; j++)
    {
        raw_output += 2 * fourier_vec_a[j] * fourier_vec_b[j];
    }

    return raw_output;
}

f_est_var Phase_based_GPR::update_estimation_essentials(const vector<Phase_based_data> &imu_data_input)
{
    // k_means_avarage
    vector<Phase_based_data> kernel_data_list;
    vector<double> phase_data;
    for (int i = 0; i < imu_data_input.size(); i++)
    {
        phase_data.push_back(imu_data_input[i].phase);
    }
    struct k_means_return k_means;
    k_means = phase_kmeans_fast(phase_data);
    for(int j = 0; j < cluster_num; j++)
    {
        struct Phase_based_data current_kernel;
        vector<double> phase_data;
        vector<double> acc_x_data;
        vector<double> acc_y_data;
        vector<double> acc_z_data;
        vector<double> agv_x_data;
        vector<double> agv_y_data;
        vector<double> agv_z_data;
        for (int i_pdata = 0; i_pdata < imu_data_input.size(); i_pdata++)
        {
            if (k_means.labels[i_pdata] == j)
            {
                phase_data.push_back(imu_data_input[i_pdata].phase);
                acc_x_data.push_back(imu_data_input[i_pdata].acc_x_imu_data);
                acc_y_data.push_back(imu_data_input[i_pdata].acc_y_imu_data);
                acc_z_data.push_back(imu_data_input[i_pdata].acc_z_imu_data);
                agv_x_data.push_back(imu_data_input[i_pdata].agv_x_imu_data);
                agv_y_data.push_back(imu_data_input[i_pdata].agv_y_imu_data);
                agv_z_data.push_back(imu_data_input[i_pdata].agv_z_imu_data);
            }
        }
        current_kernel.time_record = -1;
        current_kernel.phase = vector_mean(phase_data);
        current_kernel.acc_x_imu_data = mean_variance(acc_x_data).mean;
        current_kernel.acc_y_imu_data = mean_variance(acc_y_data).mean;
        current_kernel.acc_z_imu_data = mean_variance(acc_z_data).mean;
        current_kernel.agv_x_imu_data = mean_variance(agv_x_data).mean;
        current_kernel.agv_y_imu_data = mean_variance(agv_y_data).mean;
        current_kernel.agv_z_imu_data = mean_variance(agv_z_data).mean;
        current_kernel.sigma_acc_x = mean_variance(acc_x_data).var;
        current_kernel.sigma_acc_y = mean_variance(acc_y_data).var;
        current_kernel.sigma_acc_z = mean_variance(acc_z_data).var;
        current_kernel.sigma_agv_x = mean_variance(agv_x_data).var;
        current_kernel.sigma_agv_y = mean_variance(agv_y_data).var;
        current_kernel.sigma_agv_z = mean_variance(agv_z_data).var;
            
        kernel_data_list.push_back(current_kernel);
    }
    phase_output.clear();
    for (int i = 0; i < cluster_num; i++)
    {
        phase_output.push_back(kernel_data_list[i].phase);
    }
    ROS_WARN("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

    // update estimation essentials
    struct f_est_var est_var;
    Eigen::MatrixXd K_Gram(kernel_data_list.size(), kernel_data_list.size());
    for (int i = 0; i < kernel_data_list.size(); i++)
    {
        for (int j = 0; j < kernel_data_list.size(); j++)
        {
            K_Gram(i, j) = Fourier_kernel(kernel_data_list[i].phase, kernel_data_list[j].phase);
        }
    }
    vector<double> acc_x_list, acc_y_list, acc_z_list, agv_x_list, agv_y_list, agv_z_list;
    for (int k = 0; k < kernel_data_list.size(); k++)
    {
        acc_x_list.push_back(kernel_data_list[k].acc_x_imu_data);
        acc_y_list.push_back(kernel_data_list[k].acc_y_imu_data);
        acc_z_list.push_back(kernel_data_list[k].acc_z_imu_data);
        agv_x_list.push_back(kernel_data_list[k].agv_x_imu_data);
        agv_y_list.push_back(kernel_data_list[k].agv_y_imu_data);
        agv_z_list.push_back(kernel_data_list[k].agv_z_imu_data);
    }
    double acc_x_mean = vector_mean(acc_x_list);
    double acc_y_mean = vector_mean(acc_y_list);
    double acc_z_mean = vector_mean(acc_z_list);
    double agv_x_mean = vector_mean(agv_x_list);
    double agv_y_mean = vector_mean(agv_y_list);
    double agv_z_mean = vector_mean(agv_z_list);
    Eigen::MatrixXd Y_acc_x_osc(kernel_data_list.size(), 1), \
                    Y_acc_y_osc(kernel_data_list.size(), 1), \
                    Y_acc_z_osc(kernel_data_list.size(), 1), \
                    Y_agv_x_osc(kernel_data_list.size(), 1), \
                    Y_agv_y_osc(kernel_data_list.size(), 1), \
                    Y_agv_z_osc(kernel_data_list.size(), 1);
    for (int m = 0; m < kernel_data_list.size(); m++)
    {
        Y_acc_x_osc(m, 0) = acc_x_list[m] - acc_x_mean;
        Y_acc_y_osc(m, 0) = acc_y_list[m] - acc_y_mean;
        Y_acc_z_osc(m, 0) = acc_z_list[m] - acc_z_mean;
        Y_agv_x_osc(m, 0) = agv_x_list[m] - agv_x_mean;
        Y_agv_y_osc(m, 0) = agv_y_list[m] - agv_y_mean;
        Y_agv_z_osc(m, 0) = agv_z_list[m] - agv_z_mean;
    }
    Eigen::MatrixXd observe_noise_acc_x_diag_mat(kernel_data_list.size(), kernel_data_list.size()), observe_noise_acc_y_diag_mat(kernel_data_list.size(), kernel_data_list.size()), observe_noise_acc_z_diag_mat(kernel_data_list.size(), kernel_data_list.size());
    Eigen::MatrixXd observe_noise_agv_x_diag_mat(kernel_data_list.size(), kernel_data_list.size()), observe_noise_agv_y_diag_mat(kernel_data_list.size(), kernel_data_list.size()), observe_noise_agv_z_diag_mat(kernel_data_list.size(), kernel_data_list.size());
    for (int m = 0 ; m < kernel_data_list.size(); m++)
    {
        for (int n = 0; n < kernel_data_list.size(); n++)
        {
            if (n == m)
            {
                observe_noise_acc_x_diag_mat(m, n) = kernel_data_list[n].sigma_acc_x;
                observe_noise_acc_y_diag_mat(m, n) = kernel_data_list[n].sigma_acc_y;
                observe_noise_acc_z_diag_mat(m, n) = kernel_data_list[n].sigma_acc_z;
                observe_noise_agv_x_diag_mat(m, n) = kernel_data_list[n].sigma_agv_x;
                observe_noise_agv_y_diag_mat(m, n) = kernel_data_list[n].sigma_agv_y;
                observe_noise_agv_z_diag_mat(m, n) = kernel_data_list[n].sigma_agv_z;
            }
            else
            {
                observe_noise_acc_x_diag_mat(m, n) = 0.0;
                observe_noise_acc_y_diag_mat(m, n) = 0.0;
                observe_noise_acc_z_diag_mat(m, n) = 0.0;
                observe_noise_agv_x_diag_mat(m, n) = 0.0;
                observe_noise_agv_y_diag_mat(m, n) = 0.0;
                observe_noise_agv_z_diag_mat(m, n) = 0.0;
            }
        }
    }
    double phase_curr = imu_data_input.back().phase;
    Eigen::MatrixXd K_periodic_vec(1, kernel_data_list.size());
    for (int p = 0; p < kernel_data_list.size(); p++)
    {
        K_periodic_vec(0, p) = Fourier_kernel(phase_curr, kernel_data_list[p].phase);
    }
    Eigen::MatrixXd f_est_acc_x(1, 1);
    f_est_acc_x = K_periodic_vec * (K_Gram + observe_noise_acc_x_diag_mat).inverse() * Y_acc_x_osc;
    Eigen::MatrixXd f_est_acc_y(1, 1);
    f_est_acc_y = K_periodic_vec * (K_Gram + observe_noise_acc_y_diag_mat).inverse() * Y_acc_y_osc;
    Eigen::MatrixXd f_est_acc_z(1, 1);
    f_est_acc_z = K_periodic_vec * (K_Gram + observe_noise_acc_z_diag_mat).inverse() * Y_acc_z_osc;
    Eigen::MatrixXd f_est_agv_x(1, 1);
    f_est_agv_x = K_periodic_vec * (K_Gram + observe_noise_agv_x_diag_mat).inverse() * Y_agv_x_osc;
    Eigen::MatrixXd f_est_agv_y(1, 1);
    f_est_agv_y = K_periodic_vec * (K_Gram + observe_noise_agv_y_diag_mat).inverse() * Y_agv_y_osc;
    Eigen::MatrixXd f_est_agv_z(1, 1);
    f_est_agv_z = K_periodic_vec * (K_Gram + observe_noise_agv_z_diag_mat).inverse() * Y_agv_z_osc;
    Eigen::Vector3d acc_est;
    acc_est << f_est_acc_x(0, 0), f_est_acc_y(0, 0), f_est_acc_z(0, 0);
    Eigen::Vector3d agv_est;
    agv_est << f_est_agv_x(0, 0), f_est_agv_y(0, 0), f_est_agv_z(0, 0);

    Eigen::MatrixXd var_acc_x_est(1, 1);
    var_acc_x_est = K_periodic_vec * (K_Gram + observe_noise_acc_x_diag_mat).inverse() * K_periodic_vec.transpose();
    Eigen::MatrixXd var_acc_y_est(1, 1);
    var_acc_y_est = K_periodic_vec * (K_Gram + observe_noise_acc_y_diag_mat).inverse() * K_periodic_vec.transpose();
    Eigen::MatrixXd var_acc_z_est(1, 1);
    var_acc_z_est = K_periodic_vec * (K_Gram + observe_noise_acc_z_diag_mat).inverse() * K_periodic_vec.transpose();
    Eigen::MatrixXd var_agv_x_est(1, 1);
    var_agv_x_est = K_periodic_vec * (K_Gram + observe_noise_agv_x_diag_mat).inverse() * K_periodic_vec.transpose();
    Eigen::MatrixXd var_agv_y_est(1, 1);
    var_agv_y_est = K_periodic_vec * (K_Gram + observe_noise_agv_y_diag_mat).inverse() * K_periodic_vec.transpose();
    Eigen::MatrixXd var_agv_z_est(1, 1);
    var_agv_z_est = K_periodic_vec * (K_Gram + observe_noise_agv_z_diag_mat).inverse() * K_periodic_vec.transpose();
    Eigen::Vector3d var_acc;
    var_acc << var_acc_x_est(0, 0), var_acc_y_est(0, 0), var_acc_z_est(0, 0);
    Eigen::Vector3d var_agv;
    var_agv << var_agv_x_est(0, 0), var_agv_y_est(0, 0), var_agv_z_est(0, 0);
    Eigen::Vector3d phase_kernal;
    phase_kernal << Fourier_kernel(phase_curr, phase_curr), Fourier_kernel(phase_curr, phase_curr), Fourier_kernel(phase_curr, phase_curr);
    var_acc = -var_acc + phase_kernal;
    var_agv = -var_agv + phase_kernal;
    est_var.acc_osc = acc_est;
    est_var.agv_osc = agv_est;
    est_var.acc_var = var_acc;
    est_var.agv_var = var_agv;

    return est_var;
}

void EKF_marg::EKF_init()
{
    x.resize(7);
    x(0) = 1.0;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = 0.0;
    x(4) = 0.0;
    x(5) = 0.0;
    x(6) = 0.0;
    P = Eigen::MatrixXd::Identity(7, 7);
    Q = 1e-6 * Eigen::MatrixXd::Identity(7, 7);
    R.resize(6, 6);
    R << 10, 0, 0, 0, 0, 0,
         0, 10, 0, 0, 0, 0,
         0, 0, 10, 0, 0, 0,
         0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 1;
    R = 1e-1 * R;

    gravitational_field = Eigen::Vector3d (0, 0, 1);
    magnetic_field = Eigen::Vector3d (0.10347629, 0.99463192, 0);
}

Eigen::VectorXd EKF_marg::quaternion_multiply(const Eigen::VectorXd &q1, const Eigen::VectorXd &q2)
{
    Eigen::VectorXd q_multi(4);
    q_multi << q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3), 
               q1(0) * q2(1) + q1(1) * q2(0) + q1(2) * q2(3) - q1(3) * q2(2),
               q1(0) * q2(2) - q1(1) * q2(3) + q1(2) * q2(0) + q1(3) * q2(1),
               q1(0) * q2(3) + q1(1) * q2(2) - q1(2) * q2(1) + q1(3) * q2(0);
    return q_multi;
}

void EKF_marg::predict(const Eigen::Vector3d &gyro_data, const double &dt)
{
    // predict
    Eigen::VectorXd q(4);
    q << x(0), x(1), x(2), x(3);
    Eigen::VectorXd omega(4);
    omega << 0, gyro_data(0), gyro_data(1), gyro_data(2);
    Eigen::VectorXd q_pre_vec = q + dt * 0.5 * quaternion_multiply(q, omega);
    Eigen::Quaterniond q_pre;
    q_pre.w() = q_pre_vec(0);
    q_pre.x() = q_pre_vec(1);
    q_pre.y() = q_pre_vec(2);
    q_pre.z() = q_pre_vec(3);
    q_pre = q_pre.normalized();
    Eigen::MatrixXd F = Eigen::MatrixXd::Identity(7, 7);
    Eigen::MatrixXd F_(4, 4);
    F_ << 0, -gyro_data(0), -gyro_data(1), -gyro_data(2),
          gyro_data(0), 0, gyro_data(2), -gyro_data(1),
          gyro_data(1), -gyro_data(2), 0, gyro_data(0),
          gyro_data(2), gyro_data(1), -gyro_data(0), 0;
    Eigen::MatrixXd Iden_4 = Eigen::MatrixXd::Identity(4, 4);
    F.block<4, 4>(0, 0) = Iden_4 + 0.5 * dt * F_;
    P = F * P * F.transpose() + Q;

    x(0) = q_pre.w();
    x(1) = q_pre.x();
    x(2) = q_pre.y();
    x(3) = q_pre.z();
    x(4) = gyro_data(0);
    x(5) = gyro_data(1);
    x(6) = gyro_data(2);
}

void EKF_marg::update(const Eigen::Vector3d &acc_data, const Eigen::Vector3d &mag_data)
{
    // update
    Eigen::Quaterniond q_pre;
    q_pre.w() = x(0);
    q_pre.x() = x(1);
    q_pre.y() = x(2);
    q_pre.z() = x(3);
    Eigen::Matrix3d R_pre = q_pre.toRotationMatrix();
    Eigen::Vector3d g_pre = R_pre.transpose() * gravitational_field;
    Eigen::Vector3d m_pre = R_pre.transpose() * magnetic_field;
    Eigen::VectorXd h(6);
    h(0) = g_pre(0);
    h(1) = g_pre(1);
    h(2) = g_pre(2);
    h(3) = m_pre(0);
    h(4) = m_pre(1);
    h(5) = m_pre(2);
    Eigen::MatrixXd H(6, 7);
    H.setZero();
    H(0, 0) = 2 * (gravitational_field(0) * q_pre.w() + gravitational_field(1) * q_pre.z() - gravitational_field(2) * q_pre.y());
    H(0, 1) = 2 * (gravitational_field(0) * q_pre.x() + gravitational_field(1) * q_pre.y() + gravitational_field(2) * q_pre.z());
    H(0, 2) = 2 * (-gravitational_field(0) * q_pre.y() + gravitational_field(1) * q_pre.x() - gravitational_field(2) * q_pre.w());
    H(0, 3) = 2 * (-gravitational_field(0) * q_pre.z() + gravitational_field(1) * q_pre.w() + gravitational_field(2) * q_pre.x());

    H(1, 0) = 2 * (-gravitational_field(0) * q_pre.z() + gravitational_field(1) * q_pre.w() + gravitational_field(2) * q_pre.x());
    H(1, 1) = 2 * (gravitational_field(0) * q_pre.y() - gravitational_field(1) * q_pre.x() + gravitational_field(2) * q_pre.w());
    H(1, 2) = 2 * (gravitational_field(0) * q_pre.x() + gravitational_field(1) * q_pre.y() + gravitational_field(2) * q_pre.z());
    H(1, 3) = 2 * (-gravitational_field(0) * q_pre.w() - gravitational_field(1) * q_pre.z() + gravitational_field(2) * q_pre.y());

    H(2, 0) = 2 * (gravitational_field(0) * q_pre.y() - gravitational_field(1) * q_pre.x() + gravitational_field(2) * q_pre.w());
    H(2, 1) = 2 * (gravitational_field(0) * q_pre.z() - gravitational_field(1) * q_pre.w() - gravitational_field(2) * q_pre.x());
    H(2, 2) = 2 * (gravitational_field(0) * q_pre.w() + gravitational_field(1) * q_pre.z() - gravitational_field(2) * q_pre.y());
    H(2, 3) = 2 * (gravitational_field(0) * q_pre.x() + gravitational_field(1) * q_pre.y() + gravitational_field(2) * q_pre.z());

    H(3, 0) = 2 * (magnetic_field(0) * q_pre.w() + magnetic_field(1) * q_pre.z() - magnetic_field(2) * q_pre.y());
    H(3, 1) = 2 * (magnetic_field(0) * q_pre.x() + magnetic_field(1) * q_pre.y() + magnetic_field(2) * q_pre.z());
    H(3, 2) = 2 * (-magnetic_field(0) * q_pre.y() + magnetic_field(1) * q_pre.x() - magnetic_field(2) * q_pre.w());
    H(3, 3) = 2 * (-magnetic_field(0) * q_pre.z() + magnetic_field(1) * q_pre.w() + magnetic_field(2) * q_pre.x());

    H(4, 0) = 2 * (-magnetic_field(0) * q_pre.z() + magnetic_field(1) * q_pre.w() + magnetic_field(2) * q_pre.x());
    H(4, 1) = 2 * (magnetic_field(0) * q_pre.y() - magnetic_field(1) * q_pre.x() + magnetic_field(2) * q_pre.w());
    H(4, 2) = 2 * (magnetic_field(0) * q_pre.x() + magnetic_field(1) * q_pre.y() + magnetic_field(2) * q_pre.z());
    H(4, 3) = 2 * (-magnetic_field(0) * q_pre.w() - magnetic_field(1) * q_pre.z() + magnetic_field(2) * q_pre.y());

    H(5, 0) = 2 * (magnetic_field(0) * q_pre.y() - magnetic_field(1) * q_pre.x() + magnetic_field(2) * q_pre.w());
    H(5, 1) = 2 * (magnetic_field(0) * q_pre.z() - magnetic_field(1) * q_pre.w() - magnetic_field(2) * q_pre.x());
    H(5, 2) = 2 * (magnetic_field(0) * q_pre.w() + magnetic_field(1) * q_pre.z() - magnetic_field(2) * q_pre.y());
    H(5, 3) = 2 * (magnetic_field(0) * q_pre.x() + magnetic_field(1) * q_pre.y() + magnetic_field(2) * q_pre.z());

    Eigen::MatrixXd S = H * P * H.transpose() + R;
    Eigen::MatrixXd K = P * H.transpose() * S.inverse();

    Eigen::VectorXd z(6);
    z << acc_data(0), acc_data(1), acc_data(2), mag_data(0), mag_data(1), mag_data(2);
    x = x + K * (z - h);
    Eigen::Quaterniond q_norm;
    q_norm.w() = x(0);
    q_norm.x() = x(1);
    q_norm.y() = x(2);
    q_norm.z() = x(3);
    q_norm = q_norm.normalized();
    x(0) = q_norm.w();
    x(1) = q_norm.x();
    x(2) = q_norm.y();
    x(3) = q_norm.z();
    Eigen::MatrixXd Iden_7 = Eigen::MatrixXd::Identity(7, 7);
    P = (Iden_7 - K * H) * P;
}
