/* Generated from orogen/lib/orogen/templates/tasks/Task.cpp */

#include <boost/random/normal_distribution.hpp>
#include "Task.hpp"

using namespace imumodel;

Task::Task(std::string const& name, TaskCore::TaskState initial_state)
    : TaskBase(name, initial_state)
{
}

Task::Task(std::string const& name, RTT::ExecutionEngine* engine, TaskCore::TaskState initial_state)
    : TaskBase(name, engine, initial_state)
{
}

Task::~Task()
{
}

/// The following lines are template definitions for the various state machine
// hooks defined by Orocos::RTT. See Task.hpp for more detailed
// documentation about them.

bool Task::configureHook()
{
    if (! TaskBase::configureHook())
        return false;
    
    /** get configuration from ports */
    Configuration config;
    config.seed = _seed.value();
    config.dt = _dt.value();
    config.Dacc = _Dacc.value();
    config.Dgyro = _Dgyro.value();
    config.accrw = _accrw.value();
    config.accrrw = _accrrw.value();
    config.accbias = _accbias.value();
    config.abeta = _abeta.value();
    config.gyrorw = _gyrorw.value();
    config.gyrorrw = _gyrorrw.value();
    config.gyrobias = _gyrobias.value();
    config.gbeta = _gbeta.value();

    model.setConfiguration( config );
    model.init();

    /** Check if there is connection in the input port **/
    if (!_imuin.connected())
    {
       std::cout << "configure:: Any port connected" <<"\n";
       inport_connected = false;
    }
    else
    {
       inport_connected = true;
    }

    /** Ask if the Task is periodic and is the period is the same than the dt config value **/
    double dt = model.getConfiguration().dt;
    std::cout << "configure:: bandwidth(" << dt <<")\n";


    return true;
}
bool Task::startHook()
{
    if (! TaskBase::startHook())
        return false;

    model.reset();

    std::cout << "start:: Model initialized" <<"\n";

    return true;
}
void Task::updateHook()
{
    base::samples::IMUSensors imu_data; /**< IMU input in the port */
    base::Time timestamp = base::Time::now();

    /** Read the coming values from the inport **/
    if (inport_connected)
	_imuin.readNewest(imu_data);
    else
    {
	imu_data.acc[0] = 0.0;
	imu_data.acc[1] = 0.0;
	imu_data.acc[2] = 9.81;
	imu_data.gyro = Eigen::Matrix<double,ImuError::NUMAXIS,1>::Zero();
    }

    model.step();
    model.addNoise( imu_data );

    _imuout.write(imu_data); // write in the port
}
void Task::errorHook()
{
    TaskBase::errorHook();
    std::cout << "errorHook()" <<"\n";
}
void Task::stopHook()
{
    TaskBase::stopHook();
    std::cout << "stopHook()" <<"\n";
}
void Task::cleanupHook()
{
    TaskBase::cleanupHook();
    std::cout << "cleanupHook()" <<"\n";
}
