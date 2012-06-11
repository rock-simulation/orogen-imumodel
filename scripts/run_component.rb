#! /usr/bin/env ruby

require 'orocos'
require 'orocos/log'
require 'vizkit'

include Orocos

#Initializes the CORBA communication layer
Orocos.initialize

Orocos.run 'imumodel' do |p|
   
   # log all the output ports
   Orocos.log_all_ports 
    
   Orocos.conf.load_dir('/home/likewise-open/DFKI/jhidalgocarrio/iMoby/iMoby-dev/simulation/orogen/imumodel/configuration')
   imumodel = p.task 'imumodel'
   Orocos.conf.apply(imumodel, ['default', 'imt30'], :override => true)

   imumodel.configure
   imumodel.start
   while true
      sleep (0.1)
   end
   
end