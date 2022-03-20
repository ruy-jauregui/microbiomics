module Argu

  attr_reader :lis
  
  def initialize
   @lis = Array.new 
  end
  
	def eval_stdin
		if ARGV.length > 0
		allarg = ARGV.join(" ")
		if ARGV[0] == '-l'
			@lis = `ls`.split(/\n/)
		elsif ARGV[0] == '-f'
			@lis = IO.readlines(ARGV[1])
			@lis.each{|l| l.chomp!}
		elsif ARGV[0] =~ /\*/
			@lis = `ls #{allarg}`.split(/\n/)
		else
			@lis = ARGV
		end
		else
			puts "No args given, exiting...\n"
			exit
		end
	end

	def hal(ha,ke,va) ##a hash with arrays as values
		if ha[ke]
		  ha[ke].push(va)
		else
		  ha[ke] = Array.new
		  ha[ke].push(va)
		end
	end

end
