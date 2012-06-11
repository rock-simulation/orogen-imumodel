#! /usr/bin/env ruby

require 'orocos'
require 'orocos/log'
require 'vizkit'
require 'plotData'

include Orocos

if ARGV.size < 1 then 
    puts "usage: process_logs.rb ru"
    exit
end

def positionRbs ( data ) 
  position = [data.position.x, data.position.y, data.position.z	] 
  position
end 

def errorRbs ( data ) 
  error = [2*Math.sqrt(data.cov_position.data[0]),2*Math.sqrt(data.cov_position.data[4])]
  error
end 

def standartDeviation ( data ) 
  standart_deviation = [Math.sqrt(data.data[0]),Math.sqrt(data.data[4]),Math.sqrt(data.data[8])]
  standart_deviation
end 

def register3Axis ( title, x_axis, y_axis ) 
    plot = DataPlot.new()	
    plot.register2D( :x, {:title => "roll", :lt =>"l lt 1"} )
    plot.register2D( :y, {:title => "pitch", :lt =>"l lt 2"} )
    plot.register2D( :z, {:title => "yaw", :lt =>"l lt 3"} )
    plot.setTitle(title, "Helvetica,14")
    plot.setXLabel(x_axis, "Helvetica,14")
    plot.setYLabel(y_axis, "Helvetica,14")
    plot
end 

def dt?( data ) 
    if @init_time == 0.0 
	@init_time = data.time.to_f
    end
    dt = data.time.to_f - @init_time
    dt
end


def register3ExtraAxis(plot, name) 
    plot.register2D( :a, {:title => "roll #{name}", :lt =>"l lt 4"} )
    plot.register2D( :b, {:title => "pitch #{name}", :lt =>"l lt 5"} )
    plot.register2D( :c, {:title => "yaw #{name}", :lt =>"l lt 6"} )
end

def register3ExtraAxisTwo(plot, name) 
    plot.register2D( :u, {:title => "roll #{name}", :lt =>"l lt 7"} )
    plot.register2D( :v, {:title => "pitch #{name}", :lt =>"l lt 8"} )
    plot.register2D( :w, {:title => "yaw #{name}", :lt =>"l lt 9"} )
end

def plotArrow(plot, data)
    quaternion = Quaternion.new(data.orientation.w, data.orientation.x, data.orientation.y, data.orientation.z)
    pos = positionRbs( data )
    dcm = quaternion.q_to_dcm
    plot.arrow(pos,[pos[0]+dcm[0,0] *0.2,pos[1]+dcm[1,0] *0.2, pos[2]]) 
end

BASE_DIR = File.expand_path('..', File.dirname(__FILE__))
ENV['PKG_CONFIG_PATH'] = "#{BASE_DIR}/build:#{ENV['PKG_CONFIG_PATH']}"

Orocos.initialize

    # log all the output ports
    Orocos.log_all_ports 

    # get the invidual tasks
    log_replay = Orocos::Log::Replay.open( ARGV[0]) # This log_replay is for the logs in the folder. For the Xsens velocities and acceleration
			
    # Plot for the imu velocities (Xsens acc)
    # Plot for the orientation (Xsens versus Orientation estimator)
    plot_acc = register3Axis("Asguard acceleration", "Time (s)", "Acceleration (m/s^2)")
    plot_key = Array.new
    plot_key << :x
    plot_key << :y
    plot_key << :z
    
    sample = 0 
   
    dt = 0.0 
    init_time = 0.0
    max_accx = 0.0
    max_accy = 0.0
    max_accz = 0.0
    
    
    #log of the Xsens imu orientation
    log_replay.imumodel.imuout.connect_to :type => :buffer,:size => 2000 do|data,name|
      
      if init_time == 0.0 
	  init_time = data.time.to_f
	  max_accx = (data.acc[0]).abs
	  max_accy = (data.acc[1]).abs
	  max_accz = (data.acc[2]).abs
      end
      
#       p (data.acc[0]).abs
      
      dt = data.time.to_f - init_time
      
      if max_accx < (data.acc[0]).abs
	max_accx = (data.acc[0]).abs
      end
      if max_accy<data.acc[1].abs
	max_accy=data.acc[1].abs
      end
      if max_accz<data.acc[2].abs
	max_accz=data.acc[2].abs
      end
      
      plot_acc.addData(  :x, [dt, data.acc[0]] ) if data
      plot_acc.addData(  :y, [dt, data.acc[1]] ) if data
      plot_acc.addData(  :z, [dt, data.acc[2]] ) if data
      
      data
    
    end
    
    log_replay.align( :use_sample_time )
    control = Vizkit.control log_replay

    Vizkit.exec
  
    plot_acc.show()
    
    p max_accx
    p max_accy
    p max_accz