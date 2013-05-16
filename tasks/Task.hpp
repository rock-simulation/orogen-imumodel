/* Generated from orogen/lib/orogen/templates/tasks/Task.hpp */

#ifndef IMUMODEL_TASK_TASK_HPP
#define IMUMODEL_TASK_TASK_HPP

#include "imumodel/TaskBase.hpp"
#include <boost/random.hpp>

namespace imumodel {

    /*! \class Task 
     * \brief The task context provides and requires services. It uses an ExecutionEngine to perform its functions.
     * Essential interfaces are operations, data flow ports and properties. These interfaces have been defined using the oroGen specification.
     * In order to modify the interfaces you should (re)use oroGen and rely on the associated workflow.
     * 
     * \details
     * The name of a TaskContext is primarily defined via:
     \verbatim
     deployment 'deployment_name'
         task('custom_task_name','imumodel::Task')
     end
     \endverbatim
     *  It can be dynamically adapted when the deployment is called with a prefix argument. 
     */
    
    /** General defines **/
    #ifndef OK
    #define OK	0  /**< Integer value in order to return when everything is all right. */
    #endif
    #ifndef ERROR
    #define ERROR -1  /**< Integer value in order to return when an error occured. */
    #endif
    
    /** Sensors constant parameters **/
    #ifndef NUMAXIS
    #define NUMAXIS 3 /**< Number of axis sensed by the sensors **/
    #endif
    
    class Task : public TaskBase
    {
	friend class TaskBase;
	
    protected:
       bool inport_connected;
       base::Vector2d xax; /**< State vector for acc x axis model*/
       base::Vector2d xay; /**< State vector for acc y axis model*/
       base::Vector2d xaz; /**< State vector for acc z axis model*/
  
       
       base::Matrix2d Aax; /**< State Matrix for Accelerometer x axis */
       base::Vector2d Gax; /**< Input Matrix for Accelerometer x axis */
       base::Matrix2d Aay; /**< State Matrix for Accelerometer y axis */
       base::Vector2d Gay; /**< Input Matrix for Accelerometer y axis */
       base::Matrix2d Aaz; /**< State Matrix for Accelerometer z axis */
       base::Vector2d Gaz; /**< Input Matrix for Accelerometer z axis */
       base::Vector2d Ha; /**< Observation Matrix for Accelerometers */
      
       
       base::Vector2d xgx; /**< State vector for gyro x axis model*/
       base::Vector2d xgy; /**< State vector for gyro y axis model*/
       base::Vector2d xgz; /**< State vector for gyro z axis model*/

       
       base::Matrix2d Agx; /**< State Matrix for Gyro x axis */
       base::Vector2d Ggx; /**< Input Matrix for Gyro x axis */
       base::Matrix2d Agy; /**< State Matrix for Gyro y axis */
       base::Vector2d Ggy; /**< Input Matrix for Gyro y axis */
       base::Matrix2d Agz; /**< State Matrix for Gyro z axis */
       base::Vector2d Ggz; /**< Input Matrix for Gyro z axis */
       base::Vector2d Hg; /**< Observation Matrix for Gyroscopes */
       
       
       Eigen::Matrix <double,NUMAXIS,1> acc, gyros; /**< Accelerometer and gyroscopes output vectors */
       Eigen::Matrix <double,NUMAXIS,NUMAXIS> Dacc, Dgyro; /**< Deterministics error matrices for acc and gyros*/
       
      /*******************************************************************************
      * Normal distribution funcation which generates number with mean mean parameter 
      * and sigma standard deviation 
      * 
      ******************************************************************************/       
       typedef boost::mt19937 RandomGenerator;
       RandomGenerator rng; 
       double GetNormalDistri(double mean, double sigma);


    public:
        /** TaskContext constructor for Task
         * \param name Name of the task. This name needs to be unique to make it identifiable via nameservices.
         * \param initial_state The initial TaskState of the TaskContext. Default is Stopped state.
         */
        Task(std::string const& name = "imumodel::Task", TaskCore::TaskState initial_state = Stopped);

        /** TaskContext constructor for Task 
         * \param name Name of the task. This name needs to be unique to make it identifiable for nameservices. 
         * \param engine The RTT Execution engine to be used for this task, which serialises the execution of all commands, programs, state machines and incoming events for a task. 
         * \param initial_state The initial TaskState of the TaskContext. Default is Stopped state.
         */
        Task(std::string const& name, RTT::ExecutionEngine* engine, TaskCore::TaskState initial_state = Stopped);

        /** Default deconstructor of Task
         */
	~Task();

        /** This hook is called by Orocos when the state machine transitions
         * from PreOperational to Stopped. If it returns false, then the
         * component will stay in PreOperational. Otherwise, it goes into
         * Stopped.
         *
         * It is meaningful only if the #needs_configuration has been specified
         * in the task context definition with (for example):
         \verbatim
         task_context "TaskName" do
           needs_configuration
           ...
         end
         \endverbatim
         */
        bool configureHook();

        /** This hook is called by Orocos when the state machine transitions
         * from Stopped to Running. If it returns false, then the component will
         * stay in Stopped. Otherwise, it goes into Running and updateHook()
         * will be called.
         */
        bool startHook();

        /** This hook is called by Orocos when the component is in the Running
         * state, at each activity step. Here, the activity gives the "ticks"
         * when the hook should be called.
         *
         * The error(), exception() and fatal() calls, when called in this hook,
         * allow to get into the associated RunTimeError, Exception and
         * FatalError states. 
         *
         * In the first case, updateHook() is still called, and recover() allows
         * you to go back into the Running state.  In the second case, the
         * errorHook() will be called instead of updateHook(). In Exception, the
         * component is stopped and recover() needs to be called before starting
         * it again. Finally, FatalError cannot be recovered.
         */
        void updateHook();

        /** This hook is called by Orocos when the component is in the
         * RunTimeError state, at each activity step. See the discussion in
         * updateHook() about triggering options.
         *
         * Call recover() to go back in the Runtime state.
         */
        void errorHook();

        /** This hook is called by Orocos when the state machine transitions
         * from Running to Stopped after stop() has been called.
         */
        void stopHook();

        /** This hook is called by Orocos when the state machine transitions
         * from Stopped to PreOperational, requiring the call to configureHook()
         * before calling start() again.
         */
        void cleanupHook();
    };
}

#endif

