#! /usr/bin/env ruby

require 'orocos'
require 'orocos/log'
require 'vizkit'

include Orocos


if ARGV.size < 1 then 
    puts "usage: run_component.rb <data_log_directory>"
    exit
end

#Initializes the CORBA communication layer
Orocos.initialize

Orocos.run 'imumodel::Task'=> 'imumodel' do

    # log all the output ports
    Orocos.log_all_ports
    Orocos.conf.load_dir('../configuration/')

    imumodel = Orocos.name_service.get 'imumodel'
    Orocos.conf.apply(imumodel, ['default', 'onlyrw'], :override => true)

    # connect the tasks to the logs
    log_replay = Orocos::Log::Replay.open( ARGV[0] )

    log_replay.stim300.calibrated_sensors.connect_to(imumodel.imuin, :type => :buffer, :size => 200 )

    imumodel.configure
    imumodel.start

    # open the log replay widget
    control = Vizkit.control log_replay
    control.speed = 1 #4

    Vizkit.exec

end
