#!usr/bin/ruby --disable-gems

class Dna
  
  attr_accessor :nom, :seq, :com, :len
  
  def initialize
	@nom = ""
	@seq = ""
	@com = ""
	@len = 0
  end
  
  def setnom(aNom)
    @nom = aNom
  end
  
  def setseq(aSeq)
    @seq = aSeq
  end
  
  def setcom(aStr)
  	@com = aStr
  end
  
  def fromSl(aFile)  ##a genome sequence in a single line
   lee = File.new(aFile).gets
   lee.chomp!
   lee.downcase!
   lee.gsub!(/#/,"") ##the start char from trans-genome seq files
   if aFile =~ /\// ##if /path/included/filename
    newf = aFile.split(/\//).last
    aFile = newf
   end
   mei = aFile.gsub(/\..+$/,"")
   setnom(mei)
   setseq(lee)
   @len = lee.length
  end
  
  def fromSeq(aFile)
    mei = aFile.gsub(/\..+$/,"") ##remove extention
    mei = aFile.gsub(/-.+$/,"") ##remove positions from NC_001727-1-896.seq like
    lee = File.new(aFile).read
    lee.gsub!(/\n/,"")
    setnom(mei)
    setseq(lee)
    @len = lee.length
  end
  
  def fromGbk(af)
        myseq = ""
        ofl = false
        File.open(af).each do |l|
            if l =~ /^ORIGIN/
                  ofl = true
                  next
            end
            next unless ofl
            seql = l.chomp.gsub(/\d+|\s+/,"")
            myseq << seql if seql =~ /[atgc]/
        end
        @seq = myseq
        @len = myseq.length
        @nom = af
  end
  
  def to_faa
    lin = ">#{nom}\n#{seq}\n"
  end
  
  def print
    esc = File.new("#{@nom}","w")
    esc.print @seq
    esc.close
  end
  
  def printFaa
      esc = File.new("#{@nom.gsub(/\|/,"_")}.fa","w")
      esc.print(">#{@nom}\n")
      esc.print("#{@seq}\n")
      esc.close
  end
  
  def perm  ##random permutation of seq
    newseq = @seq.scan(/./).shuffle.join
    @seq = newseq
  end
  
  def ran(aLen) ##generates random sequence of aLen (integer) length
    sou = ["a","t","g","c"]
    rseq = ""
    for j in 0..aLen
      rseq << sou[rand(4)]
    end
    @nom = "random"
    @seq = rseq
  end
  
  def mask
    @seq.tr!('rymkdhvbswn', 'aaaaaaaaaaa')
  end
  
  def to_curva  ##format for running the curva program
    #mask
	 @seq.gsub!(/\n/,"")
    j = @seq.length
    newseq = ""
    z = 0
    until z > j 
      newseq << "#{@seq[z..(z+59)]}\n"
      z += 60
    end
    @seq = newseq
  end
  
  def nucc ###Nucleotide composition
    nuc = "a t g c".split
    sal = ""
    di = Hash.new
    di.default = 0
    sl = @seq.length
    sn = @seq.downcase.split
    sn.each{|k|
      next if k =~ /[rmykwbdhvsn]/
      di[k] = di[k] + 1
    }
    to = sl.to_f
    nuc.each{|n|
     if di[n]
       fre = (di[n] / to)
       sal << "#{n}\t#{fre}\n"
     else
       sal << "#{n}\t0.0\n"
     end
    }
    sal
  end
  
  def dinuc  ##dinucleotide comp.
    nuc = "aa ac at ag ca cc ct cg ta tc tt tg ga gc gt gg".split
    sal = ""
    dih = Hash.new(0)
    sl = @seq.length
    for i in (0..(sl - 1))
      din = @seq[i..(i+1)]
      next if din =~ /[rmykwbdhvsn]/ 
      if din =~ /[atgc]/
            dih[din]  += 1 if din.length  == 2
      end
    end
    to = sl.to_f
    nuc.each{|n|
     if di[n]
       fre = di[n] / to
       sal << "#{n}\t#{fre}\n"
     else
       sal << "#{n}\t0.0\n"
     end
    }
    sal
  end
  
  def trinuc
      a = %w(a t g c)
      tria = [] ##all trinucleotides.
      a.repeated_permutation(3).each{|e| tria.push(e.join)}
      trih = Hash.new(0)
      sl = @seq.length
      for i in (0..(sl - 2))
            trn = @seq[i..(i+2)] #get a trn
            next if trn =~ /[rmykwbdhvsn]/ #supposing it is unmasked, skip the offending trinucleotide
            trih[trn] += 1 if trn.length == 3
      end
      sal = ""
      tria.each{|k|
            if trih[k]
                 fre = trih[k] / (sl - 2)
                 #puts "#{k},#{fre}"
                 sal << "#{k},#{fre}\n"
            else
                  #puts "#{k},#{0}"
                  sal << "#{k},0\n"
            end
            
      }
      sal
  end
  
  def gc ##evaluate gc %
      len = @seq.length
      nel = len.to_f
      gt = @seq.downcase.count "gc"
      per = (gt / nel) * 100
	per
  end
  
  def slidegc(step) ##evaluate gc% in a sliding window
  	sal = String.new
  	fin = @seq.length - step - 1
	1.upto(fin) {|j|
		frag = @seq[j..(j+step)]
		num = frag.count "gc"
		gc = (num/step) * 100
		sal << "#{gc}\n"
	}
	sal
  end
  
	def stepgc(step) 
  		sal = String.new
  		fin = @seq.length - step - 1
		j = 1
		while (j <= fin)
			frag = @seq[j..(j+step)]
			frag.downcase!
			#puts "#{frag}\n\n"
			num = frag.count "gc"
			#puts "num is #{num}\n\n"
			gc = (num.to_f / step.to_f) * 100
			sal << "#{gc}\n"
			j += step
		end
		sal
	end
	
	def stepslidegc(step,win)
		sal = String.new
  		fin = @seq.length - win - 1
		j = 1
		while (j <= fin)
			frag = @seq[j..(j+win)]
			frag.downcase!
			num = frag.count "gc"
			gc = (num.to_f / win.to_f) * 100
			sal << "#{gc}\n"
			j += step
		end
		sal
	end

end
class String
      def revcomp #to make reverse complement of @seq
            self.tr('aAtTgGcCrRyYmMkKdDhHvVbBsSwWnN', 'tTaAcCgGyYrRkKmMhHdDbBvVsSwWnN').reverse
      end   
end
