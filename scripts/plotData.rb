require 'gnuplot'

class DataPlot


    def initialize(  )
	@plots = Array.new
	@objects = Array.new
	@sets = Array.new 
	@figure = Array.new 
	@dimension = 0 
	@hash_order = Hash.new
    end

    def registered?(sym)
        @plots.has_key?(sym)
    end
    
    def id( sym ) 
	@hash_order[sym]
    end 
    
    def registery( sym, params )
	@hash_order[sym] = @hash_order.size
	
	@plots[id(sym)] = Hash.new
	@plots[id(sym)].merge!( params )
    end
 
    def register1D( sym, params )
	@dimension = 1
	registery( sym, params )
	@plots[id(sym)][:data] = [[]] 
    end
    
    def register2D( sym, params )
	@dimension = 2
	registery( sym, params )
	@plots[id(sym)][:data] = [ [], [] ] 
    end
    
    def register3D( sym, params )
	@dimension = 3
	registery( sym, params )
	@plots[id(sym)][:data] = [ [], [], [] ] 
    end
    
    #type exampel "linespoints" or "points" 
    def plotType(sym, type)
	@plots[id(sym)][:lt] = type
    end 
 
    
    def addData( sym, data)
	if @dimension == 1
	    @plots[id(sym)][:data][0] << data
	elsif @dimension >= 2
	    @plots[id(sym)][:data][0] << data[0]
	    @plots[id(sym)][:data][1] << data[1]
	end
	if @dimension == 3
	    @plots[id(sym)][:data][2] << data[2]
	end 

    end
    
    def arrow(from, to, options = nil)
	
	if !options 
	    option = ""
	end
	
	obj_idx = @objects.size + 1
	
	if @dimension == 3
	    @objects << "arrow #{obj_idx} from #{from[0]}, #{from[1]}, #{from[2]} to #{to[0]}, #{to[1]}, #{to[2]} #{options}"
	else 
	    @objects << "arrow #{obj_idx} from #{from[0]}, #{from[1]} to #{to[0]}, #{to[1]} #{options}"
	end
	
    end 
    
    def ellipse(center, size, options = nil)
      
	if !options
	    options = "" 
	end
	
	obj_idx = @objects.size + 1
	
	@objects << "object #{obj_idx} ellipse center #{center[0]}, #{center[1]} size #{size[0]},#{size[1]} #{options}" 
	
    end 
    
    def clearObjects()
	@objects.clear
    end 
    
    def generalPurposeSet( set ) 
	@sets << set 
    end 
    
    def set( set ) 
	@sets << set 
    end  
    
    #Font example - Helvetica,14
    def setTitle ( title, font = nil ) 
      	if !font 
	    font = "Helvetica,14" 
	end
	@sets << "title '#{title}' font '#{font}'"
    end 
    
    def setGrid ()
	@sets << "grid" 
    end
    
    def setFigure( name, terminal_size) 
	@figure.clear
	@figure << "terminal png large size #{terminal_size[0]},#{terminal_size[1]}"
	@figure << "output '#{name}'"
    end 
    
    #Font example - Helvetica,14
    def setXLabel( name, font = nil )
	if !font 
	    font = "Helvetica,14" 
	end
      	@sets << "xlabel '#{name}' font '#{font}'"
    end 
    
    def setYLabel( name, font = nil )
	if !font 
	    font = "Helvetica,14" 
	end
	@sets << "ylabel '#{name}' font '#{font}'"
    end 
    
    def setZLabel( name, font = nil)
	if !font 
	    font = "Helvetica,14" 
	end
      	@sets << "zlabel '#{name}' font '#{font}'"
    end 
  
    #right bottom  Helvetica,14   1.5
    def setLabel ( position, font, spacing ) 
	@sets <<  "key spacing #{spacing}"
	@sets <<  "key font '#{font}'"
	@sets <<  "key #{position}"
    end 
    
    def setXRange ( range )
        @sets << "xrange [#{range[0]}:#{range[1]}]"
    end 

    def setYRange ( range )
        @sets << "yrange [#{range[0]}:#{range[1]}]"
    end 
    
    def setZRange ( range )
        @sets << "zrange [#{range[0]}:#{range[1]}]"
    end 

    def save( file_name )
	File.open( file_name, 'w') do |file|
	    gnuplots( file )
	end
    end
    
    def show( )
	Gnuplot.open do |gp|
	    gnuplots( gp )
	end
    end

      def gnuplots( gp )
	    if @dimension == 2 || @dimension == 1 
		Gnuplot::Plot.new( gp ) do |plot|
		    gnuplot ( plot ) 
		end
	    elsif @dimension == 3
		#Gnuplot::SPlot.new( gp ) do |plot|
		Gnuplot::Plot.new( gp ) do |plot|
		    gnuplot ( plot ) 
		end
	    end
    end
	       
    def gnuplot ( plot ) 
	
	plots = Array.new
	
	for i in 0..@sets.size-1
	    plot.set @sets[i]
	end
	for i in 0..@objects.size-1
	    plot.set @objects[i]
	end
	for i in 0..@figure.size-1
	    plot.set @figure[i]
	end
	for i in 0..@plots.size-1
	    if @plots[i][:data][0].length > 0
		
		plots << Gnuplot::DataSet.new( @plots[i][:data] ) do |ds|
		    @plots[i].each do |key, value|
			with = "" 
			if @plots[i][:with]
			    with = with + " #{@plots[i][:with]} "
			end
			if @plots[i][:lt]
			    with = with + " #{@plots[i][:lt]} "
			end
			if with != ""
			    ds.with = with
			end
			
			if key != :data && key !=:lt && key !=:with 
			    ds.send("#{key}=", value)
			end
		    end
		end
	    end
	    
	end
	plot.data = plots
    end 
	       
end
