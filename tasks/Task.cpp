/* Generated from orogen/lib/orogen/templates/tasks/Task.cpp */

#include <boost/random.hpp>
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
   base::Vector3d abeta = _abeta.get();
   base::Vector3d accrrw = _accrrw.get();
   base::Vector3d accbias = _accbias.get();
   base::Vector3d gbeta = _gbeta.get();
   base::Vector3d gyrorrw = _gyrorrw.get();
   base::Vector3d gyrobias = _gyrobias.get();
   double dt = _dt.value();
    
    if (! TaskBase::configureHook())
        return false;
    
    
    
    /** Fill the matrices with the configure values for the accelerometer model X axis **/
     Aax(0,0) = 1.0;
     Aax(0,1) = 0.0;
     Aax(1,0) = dt;
     Aax(1,1) = 1-(abeta[0]*dt);
     
     Gax (0,0) = accrrw[0]* abeta[0];
     Gax (1,0) = (accrrw[0]*sqrt(dt/3))+ abeta[0]*accbias[0];
     
     std::cout << "Xaxis\n";
     std::cout << Aax<<"\n";
     std::cout << Gax<<"\n";
     
    
    /** Fill the matrices with the configure values for the accelerometer model Y axis **/
    Aay(0,0) = 1.0;
    Aay(0,1) = 0.0;
    Aay(1,0) = dt;
    Aay(1,1) = 1-(_abeta.value()(1)*_dt.value());
    
    Gay (0,0) = accrrw[1]* _abeta.value()(1);
    Gay (1,0) = (accrrw[1]*sqrt(dt/3))+ abeta[1]*accbias[1];
    
    std::cout << "Yaxis\n";
     std::cout << Aay<<"\n";
     std::cout << Gay<<"\n";
    
    /** Fill the matrices with the configure values for the accelerometer model Z axis **/
    Aaz(0,0) = 1.0;
    Aaz(0,1) = 0.0;
    Aaz(1,0) = dt;
    Aaz(1,1) = 1-(_abeta.value()(2)*_dt.value());
    
    Gaz (0,0) = accrrw[2]* _abeta.value()(1);
    Gaz (1,0) = (accrrw[2]*sqrt(dt/3))+ abeta[2]*accbias[2];
    
    std::cout << "Zaxis\n";
     std::cout << Aaz<<"\n";
     std::cout << Gaz<<"\n";
    
    /** Measurement matrix **/
    Ha << 0,1;
    
    
    /** Fill the matrices with the configure values for the gyroscope model X axis **/
    Agx(0,0) = 1.0;
    Agx(0,1) = 0.0;
    Agx(1,0) = dt;
    Agx(1,1) = 1.0-(gbeta[0]*dt);
    
    Ggx (0,0) = gyrorrw[0]* gbeta[0];
    Ggx (1,0) = (gyrorrw[0]*sqrt(dt/3)) + gbeta[0]*gyrobias[0];
    
    std::cout << "Xaxis\n";
     std::cout << Agx<<"\n";
     std::cout << Ggx<<"\n";
    
    /** Fill the matrices with the configure values for the gyroscope model Y axis **/
    Agy(0,0) = 1.0;
    Agy(0,1) = 0.0;
    Agy(1,0) = dt;
    Agy(1,1) = 1-(gbeta[1]*dt);
    
    Ggy (0,0) = gyrorrw[1]* gbeta[1];
    Ggy (1,0) = (gyrorrw[1]*sqrt(dt/3)) + gbeta[1]*gyrobias[1];
    
    std::cout << "Yaxis\n";
     std::cout << Agy<<"\n";
     std::cout << Ggy<<"\n";
    
    /** Fill the matrices with the configure values for the gyroscope model Z axis **/
    Agz(0,0) = 1.0;
    Agz(0,1) = 0.0;
    Agz(1,0) = dt;
    Agz(1,1) = 1-(gbeta[2]*dt);
    
    Ggz (0,0) = gyrorrw[2]* gbeta[1];
    Ggz (1,0) = (gyrorrw[2]*sqrt(dt/3)) + gbeta[2]*gyrobias[2];
    
    std::cout << "Zaxis\n";
     std::cout << Agz<<"\n";
     std::cout << Ggz<<"\n";
     
    /** Measurement matrix **/
    Hg << 0,1;
    
    /** deterministic matrices to zero **/
    Dacc = Eigen::Matrix<double,NUMAXIS,NUMAXIS>::Zero();
    Dgyro = Eigen::Matrix<double,NUMAXIS,NUMAXIS>::Zero();
    
    
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
    if ((TaskContext::getPeriod() != dt) && (TaskContext::getPeriod() != 0.0))
    {
	 std::cout<< "configure(warning): the task has two periods values\n";
	 TaskContext::setPeriod(dt);
	 std::cout<< "configure: the task is configured to("<<dt<<") seconds period\n";
    }
      else
	 std::cout << "configure:: period(" << TaskContext::getPeriod() <<")\n";

    
    return true;
}
bool Task::startHook()
{
    if (! TaskBase::startHook())
        return false;
    
    
    /** Reset the state vectors and the output vectors **/
    xax << 0.0,0.0;
    xay << 0.0,0.0;
    xaz << 0.0,0.0;
    
    xgx << 0.0,0.0;
    xgy << 0.0,0.0;
    xgz << 0.0,0.0;
    
    acc = Eigen::Matrix<double,NUMAXIS,1>::Zero();
    gyros = Eigen::Matrix<double,NUMAXIS,1>::Zero();
    
    Dacc = _Dacc.value();
    Dgyro = _Dgyro.value();
    
    
   std::cout << "start:: Model initialized" <<"\n";

    
    return true;
}
void Task::updateHook()
{
   base::samples::IMUSensors imu_in; /**< IMU input in the port */
   base::samples::IMUSensors imu_out; /**< IMU final output to the port */  
   base::Time timestamp = base::Time::now();
       
   //TaskBase::updateHook();
//    std::cout<< "Time: " <<timestamp.toString() <<"\n";
   
   /** Read the coming values from the inport **/
   if (inport_connected)
      _imuin.readNewest(imu_in);
   else
   {
      imu_in.acc[0] = 0.0;
      imu_in.acc[1] = 0.0;
      imu_in.acc[2] = 9.81;
      imu_in.gyro = Eigen::Matrix<double,NUMAXIS,1>::Zero();
   }
    
    /** White noise for the Velocty Random Walk for accelerometers **/
    base::Vector3d stdacc;
    stdacc[0] = _accrw.value()[0]/sqrt(_dt.value());
    stdacc[1] = _accrw.value()[1]/sqrt(_dt.value());
    stdacc[2] = _accrw.value()[2]/sqrt(_dt.value());
    std::cout << "Standard deviation(acc)\n" <<stdacc <<"\n";
  
    /** Perform the model step for Accelerometers Xaxis **/
    xax = Aax * xax + Gax * GetNormalDistri(0,1.0);
    acc[0] = GetNormalDistri(0,stdacc[0]) + (Ha.transpose() * xax); // White noise + stochastic noise
    
    /** Perform the model step for Accelerometers Yaxis **/
    xay = Aay * xay + Gay * GetNormalDistri(0,1.0);
    acc[1] = GetNormalDistri(0,stdacc[1]) + (Ha.transpose() * xay); // White noise + stochastic noise
    
    /** Perform the model step for Accelerometers Zaxis **/
    xaz = Aaz * xaz + Gaz * GetNormalDistri(0,1.0);
    acc[2] = GetNormalDistri(0,stdacc[2]) + (Ha.transpose() * xaz); // White noise + stochastic noise
    
//     std::cout << "xax\n" << xax <<"\n";
//     std::cout << "xay\n" << xay <<"\n";
//     std::cout << "xaz\n" << xaz <<"\n";
//     std::cout << "ACC\n" << acc <<"\n";
    
    /** White noise for the Angular Random Walk for gyroscopes **/

    base::Vector3d stdgyro;
    stdgyro[0] = _gyrorw.value()[0]/sqrt(_dt.value());
    stdgyro[1] = _gyrorw.value()[1]/sqrt(_dt.value());
    stdgyro[2] = _gyrorw.value()[2]/sqrt(_dt.value());
    std::cout << "Standard deviation(gyro)\n" <<stdgyro <<"\n";
    
    /** Perform the model step for gyroscopes Xaxis **/
    xgx = Agx * xgx + Ggx * GetNormalDistri(0,1.0);
    gyros[0] = GetNormalDistri(0, stdgyro[0]) + (Hg.transpose() * xgx); // White noise + stochastic noise
    
    /** Perform the model step for gyroscopes Yaxis **/
    xgy = Agy * xgy + Ggy * GetNormalDistri(0,1.0);
    gyros[1] = GetNormalDistri(0, stdgyro[1]) + (Hg.transpose() * xgy); // White noise + stochastic noise
    
    /** Perform the model step for gyroscopes Zaxis **/
    xgz = Agz * xgz + Ggz * GetNormalDistri(0,1.0);
    gyros[2] = GetNormalDistri(0, stdgyro[2]) + (Hg.transpose() * xgz); // White noise + stochastic noise
    

//     std::cout << "xgx\n" << xgx <<"\n";
//     std::cout << "xgy\n" << xgy <<"\n";
//     std::cout << "xgz\n" << xgz <<"\n";
//     std::cout << "GYRO\n" << gyros <<"\n";
//        
    /** Include the deterministic error to the model and white it to the output ports (accelerometer and gyroscopes) **/
    imu_out.acc = (Eigen::Matrix<double,NUMAXIS,NUMAXIS>::Identity() + Dacc) * imu_in.acc + acc;
    imu_out.gyro = (Eigen::Matrix<double,NUMAXIS,NUMAXIS>::Identity() + Dgyro) * imu_in.gyro + gyros;
    imu_out.time = timestamp;
    
    _imuout.write(imu_out); // write in the port
    
//     std::cout << "acc_out\n" << imu_out.acc <<"\n";
//     std::cout << "gyro_out\n" << imu_out.gyro <<"\n";

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



double Task::GetNormalDistri(double mean, double sigma)
{
 typedef boost::normal_distribution<double> NormalDistribution;
 typedef boost::mt19937 RandomGenerator;
 typedef boost::variate_generator<RandomGenerator&,NormalDistribution> GaussianGenerator;
 
  /** Initiate Random Number generator with current time */
  static RandomGenerator rng(static_cast<unsigned> (time(0)));
 
  /* Choose Normal Distribution */
  NormalDistribution gaussian_dist(mean, sigma);
 
  /* Create a Gaussian Random Number generator
   *  by binding with previously defined
   *  normal distribution object
   */
  GaussianGenerator generator(rng, gaussian_dist);
 
  // sample from the distribution
  return generator();
}
