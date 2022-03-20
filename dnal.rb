#!/home/jaureguir/bin/ruby --disable-gems
##class to manipulate lists of sequences
require "/home/jaureguir/bin/dna.rb"
class Dnal
  
  attr_accessor :faaseq, :faaha
  
  def initialize
	@faaseq = Array.new ## DNA objects
	@faaha = Hash.new
  end
  
  def toCurva
    @faaseq.each{|v| v.to_curva}
  end
  
  def frag(pos,len) ##get only a fragment of the sequence(position, length)
  	le = len.to_i ##length
  	po = pos.to_i ##initial position
  	fp = po + le ##final position
  	@faaseq.each{|o|
  		dum = o.seq.dup
  		if fp <= dum.length
  			neu = dum[po..fp]
  			o.setseq(neu)
  		end
  	}
  end
  
  def flock(le) ##generate new objects with a fragment of the sequence, return an array of dna objects, le is the fragment length
  	neufas = Array.new
  	@faaseq.each{|o|
  		pos = 0
  		fr = 1
  		#puts "original seq #{o.seq}"
  		le = o.seq.length if le > o.seq.length
  		while pos < o.seq.length
  			break if (pos + le) > o.seq.length
  			if (pos + le) < o.seq.length
  			frag = o.seq[pos,le]
  			#puts "frag: #{frag}"
  			next if frag.length < le
  			pos += le
  			#puts "current pos #{pos}"
  			neo = Dna.new
			neo.setnom("#{o.nom}_#{fr}")
			neo.setseq(frag)
			#puts "new obj: #{p neo}"
			neufas.push(neo)
			fr += 1
			end
  		end
  	}
  	neufas
  end
  
  def toFaa(arg) #haaseq or faaha
  	lin = ""
  	if arg.class == Array
           arg.each{|o|
		lin << o.to_faa
            } 
  	elsif arg.class == Hash
            arg.each{|k,o|
		    lin << o.to_faa
            }
  	end
	lin
  end
  
  def artoFaa
      lin = ""
      @faaseq.each{|v|
            lin << v.to_faa
            lin << "\n"
      }
      lin
  end
  
  def hatoFaa
      lin = ""
      @faaha.each{|k,v|
            lin << v.to_faa      
      }
      lin
  end
  
  def fromSl(aFile)##one line seq
   sek = Dna.new
   sek.fromSl(aFile)
   @faaseq.push(sek)
  end
  
  def fromFaa(aFile,apa)
    rege = "Iwillnevereversplithteheaderofthisfastaseq" 
    if apa.length > 0
      #apa.gsub(/\|/,"\|") if apa =~ /\|/
      rege = Regexp.new(apa)
      #rege = Regexp.new(apa).escape('|') if apa =~ /\|/
    end  
    srand;
    @tm = Hash.new
    fun = Hash.new
    lee = File.open(aFile)
    he = Array.new
    j = rand(100000)
    while lee.gets
      if $_.chars.select{|i| i.valid_encoding?}.join =~ /^>/
	alli = $_.chars.select{|i| i.valid_encoding?}.join.chomp.gsub(/>/,"").split(rege)
	alli[0].gsub!(/\s+$/,"") #remove trailing spaces or tabs
	#alli[0].gsub!(/ /,"") #if this is commented, spaces are allowed in the header
	head = alli[0]
	#head = alli[3].gsub(/\s/,"_")
	fun[head] = $_.chars.select{|i| i.valid_encoding?}.join.chomp.gsub(/>/,"")
	he[0] = "#{head}::#{j}" ##to avoid adding seq to repeated headers
	j = rand(100000)
      else
	if @tm[he[0]] 
	  @tm[he[0]] << $_.chars.select{|i| i.valid_encoding?}.join.chomp ##adds a single line
	else
	  @tm[he[0]] = ""
	  @tm[he[0]] << $_.chars.select{|i| i.valid_encoding?}.join.chomp
	end
      end
    end
    #p @tm
    
    @tm.each{|x,y|
      lxa = x.split(/::/)
      sek = Dna.new
      sek.setnom(lxa[0])
      sek.setseq(y)
      sek.setcom(fun[lxa[0]])
      sek.len = y.length
      @faaseq.push(sek)
      @faaha[sek.nom] = sek
    }
  end
  
  def fromList(aFile)
    File.readlines(aFile).each{|l|
      cham = l.split ##watch out!
      sek = Dna.new
      sek.setnom(cham[0])
      sek.setseq(cham.last.chomp)
      sek.len = cham.last.length
      @faaseq.push(sek)
    }
  end
  
  def fromHash(aHash) ##in memory, id => seq
    aHash.each{|x,y|
      sek = Dna.new
      sek.setnom(x)
      sek.setseq(y)
      sek.len = y.length
      @faaseq.push(sek)
    }
  end
  
  def fromSeq(aFile) ##file.seq, in curva format
    sek = Dna.new
    sek.fromSeq(aFile)
    @faaseq.push(sek)
  end
  
  def fromGbk(af) #a *.gbk file
        sek = Dna.new
        sek.fromGbk(af)
        @faaseq.push(sek)
  end
  
  def nr ##eliminates repeated objects (same name)
  	hak = Hash.new
	tar = Array.new
	@faaseq.each{|o|
		next if hak[o.nom]
		tar.push(o)
		hak[o.nom] = 1
	}
	@faaseq = tar
  end
  
  def getsome(num) ##gets n entries at random
  	srand
  	sel = Array.new
	lar = @faaseq.length
	dumy = @faaseq.dup
	for i in (1..num)
		r = rand(lar)
		sel.push(dumy[r])
		dumy.delete_at(r)
		lar -= 1
	end
	sel
  end
  
  def doha #done by default in fromFaa, to use when fromFaa is not called!
  	@faaseq.each{|e| @faaha[e.nom] = e }
  end
  
  def reha #reset the hash keys if the name has changed
      nha = {}
      @faaha.each{|k,v| nha[faaha[k].nom] = v}
      @faaha = nha
  end
  
  def perm ##random permutation of each sequence, adds "-rand" to e.nom
  	@faaseq.each{|e| e.perm}
  end

  def headers
      return @faaha.keys
  end
  
  def modhead(opa, npa) #modifies the header of all files, old pattern, new pattern
    #the usual would be: ":|#","_"
      olp = Regexp.new(opa)
      #olp.escape('|') if olp =~ /\|/
      @faaseq.each{|d| d.nom.gsub!(olp,npa) }
  end
  
  def revcomp
      @faaseq.each{|d| rc = d.seq.revcomp; d.seq = rc}
  end
  
  def seqindex #make a hash with sequence as key and name as value
      seqindh = Hash.new
      faaseq.each{|f| seqindh[f.seq] = f.nom}
      return seqindh
  end
  
  def lensort #sort by seq length
      @faaseq.sort!{|a,b| b.len <=> a.len}
  end
  
  def remlen(an) #remove if seq shorter than an
      ai = an.to_i
      @faaseq.delete_if{|x| x.len < ai}
      lensort
  end
  
  def addhead(as) #adds a string to the beginnign of the fasta header
    faaseq.each{|f| f.nom = "#{as}_#{f.nom}"}
  end
  
  #deal with fastq files, use Dna.com to store the qual line
  def fromFastq(aq)
    fastqf = File.new(aq)
    ### THIS BELOW DEPENDS ON THE SOURCE HEADER FORMAT. ":" for illumina, "_" for PacBio
    headerpattern = Regexp.new(fastqf.readline.split(/_/)[0]) #illumina's first name in header
    fastqf.rewind #set IO stream at beginning again, line 0
    i = 0
    myd = Dna.new #will never be used, but sets myd class.
    fastqf.each do |l|
      if l=~ headerpattern #the header
				myd = Dna.new
				myd.nom = l.chomp
				i = 0
      end
      if i == 1 #the sequence
				myd.seq = l.chomp
      end
      if i == 3 #the qual line
				myd.com = l.chomp
				@faaseq.push(myd)
      end
      i += 1
    end
  end
  
  def toFastq
    lin = ""
    @faaseq.each{|f| lin << "#{f.nom}\n#{f.seq}\n+\n#{f.com}\n"} #the quality is in the com field
    return lin
  end
  
  def degap #remove gaps if the sequences were aligned.
		@faaseq.each{|f| f.seq.gsub!(/-/,"")}
	end
  
end
