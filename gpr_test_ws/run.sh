#!/bin/bash
source devel/setup.bash
roslaunch gpr gpr.launch & rosbag record -O gpr.bag /mavros/imu/data /mavros/imu/data_raw /mavros/imu/mag /mavros/global_position/raw/fix /mavros/global_position/raw/gps_vel /imu_osc /imu_no_osc
