#!/usr/bin/ruby --disable-gems
#microbiome filter method
#reads data from the pipeline results described in: Camarinha-Silva A, Jáuregui R, Chaves-Moreno D, Oxley AP, Schaumburg F, Becker K, Wos-Oxley ML, Pieper DH. Comparing the anterior nare bacterial community of two discrete human populations using Illumina amplicon sequencing. Environ Microbiol. 2014 Sep;16(9):2939-52. doi: 10.1111/1462-2920.12362. Epub 2014 Jan 16. PMID: 24354520.
require 'dnal.rb'
require 'argu.rb'

class Gut
      
      include Argu
      
      def filter(clu,fas) # unique precluster name file from mothur, precluster fasta file
            cle = clu.split(/\./) #to write the resulting files
            stamp = Time.new.strftime("%Y%m%d-%H%M%S") #time stamp for log file
            lof = File.open("#{cle[0]}-#{stamp}-log.csv","w") #the log file
            lof.print("query files: #{clu}, #{fas}\n")
            cos = readspersample(clu) # reads per sample: hash Sample name - > read No.
            oneinth = {} #will hold the one in ten thousand thresholds
            cos.each{|k,v| oneinth[k] = v/10000} #the one in ten thousand thresholds
            grandtot = 0 #total number of reads in all samples
            cos.values.each{|v| grandtot += v }
            grandthr = 0.01 #the threshold, in %, on the abundance of the overall community
            puts "total read No: #{grandtot}"
            mincom = (grandtot.to_f * grandthr) / 100.0 #minimal abundance in the whole set, including all samples
            #or override the overall community filter here, if you know how much you want:
            mincom = 1#40 #at least 40 reads in all, 1: no filter
            puts "#{grandthr}% of overall community: filter is #{mincom.to_i}"
            lof.print("#{grandthr}% of overall community: filter is #{mincom.to_i}\n")
            ####
            #the thr1, testing if the rep is more than 1% of a single sample, is tested below, see variable fla1.
            ####
            thrh = {} #threshold on abundance in each sample, to combine with the No of samples that have it
            samabper = 0.2
            cos.each{|k,v| thrh[k] = (v.to_f * samabper) / 100.0 ; puts "#{k} tot:#{v} #{samabper}% abundance thr: #{thrh[k]}"
                  lof.print("#{k}  tot:#{v} #{samabper}% abundance thr: #{thrh[k]}\n")
            }#abundance per sample, to be used in combination with the % of samples that have it above this thr.
            
            #thr5 is the No. of samples that have a given read, if the read is in more than thr5, it is taken without abundance check.
            samper = 5.0 #percentage of samples that have a given rep,  -> no abundance check
            
            thr5 = (cos.keys.length.to_f * samper) / 100.0 #total sample number based on #{samper} %, beware! here that % of samples must have the read for the read to be kept
            thr5 = 2.0 if thr5 < 2.0 #floor threshold, the read is in at least 5 samples, -> no abundance check
            puts "OTUs will be kept if their local abundance is >#{grandthr * 100}% (in any given sample), and following these thresholds:"
            puts "Sample number threshold, #{samper} %, if an OTU is in #{thr5} samples, it is kept with no abundance check"
            lof.print("Sample number threshold, #{samper} %, if an OTU is in #{thr5} samples with any abundance, it is kept \n")
            #thr2 is the minimal sample No with read abundance bigger than thrh[k]
            sammanyper = 2.0 # this value in x% of samples must have the rep with abundance higher than #{samabper} %            
            thr2 = (cos.keys.length.to_f * sammanyper ) / 100.0 # n% of samples, to implement in combination with the abundance threshold (thrh)
            thr2 = 2.0 if thr2 < 2.0 #floor threshold
            puts "the thr2, or #{sammanyper}% of total samples, each with with >= #{samabper}% read abundance in each sample is #{thr2}"
            lof.print("the thr2, or #{sammanyper} % of total samples, each with with >= #{samabper}% read abundance in each sample is #{thr2}\n")
            maxthr = 1 #minimal copy No. threshold. At least one sample must have maxthr reads
            lof.print("minimal copy No per sample is set to #{maxthr}\n")
            repcps = {} #key is a representative seq id, value is an array of sample -> total count (adding up the clustered copies)
            keep = [] #will hold seq id's to keep, having passed ALL FILTERS<---!! (wow!)
            lin = "Index,seq,Total_copy_No,in_Sample_No,#{cos.keys.sort.join(",")}\n"
            lcn = "Index,seq,Total_copy_No,in_Sample_No,#{cos.keys.sort.join(",")}\n"
            faline = "" #the fasta line
            myd = Dnal.new
            myd.fromFaa(fas," ") ##<---------Warinig: BEWARE OF THE HEADER PATTERN!!! Here it splits by " " and keeps only header[0]
            index = 1 #this will be the sequence name now
            File.open(clu).each do |s| #precluster name line
                  (he, re) = s.chomp.split(/\s+/) #representative seq, list of clustered seq ids
                  tcc = re.split(/,/).length.to_f #phylotype total copy count in all samples, add all copies from different clustered seqs
                  tsc = 0 #total sample count,
                  cps = {} #copies per sample
                  alls = [] #all samples pool together here
                  re.split(/,/).each{|sa| alls.push(sa.split(/\||_/)[0])} #push sample names into alls
                  alls.uniq.each{|e| cps[e] = alls.count(e).to_f } #sample name -> No of times it appears in the precluster line
                  tsc = alls.uniq.length.to_f #total sample count, total number of samples that have this read, alls has sample names
                  fla5 = false
                  minreadnum = false #minimal read number, will turn true if one sample has >= than maxthr reads
                  
                  fla1 = false #will turn true if one sample has the read in abundance >= 1%
                  lowabc = 0.0 #low abundance count (samples with more than 0.1%)
                  countthr5 = 0 #count samples that pass thr5 abundance
                  lowfla = false
                  countone = 0 #count of samples with read abundance >1%
                  countzero = 0 #count of samples with abundance > 0.1%
                  cps.each do |k,v| #sampe name -> read number in this line, cps has copies per sample
                        minreadnum = true if v >= maxthr #keep if at least maxthr reads are in one sample
                        if cos[k] > 9999.0 #cos[k] is total reads in sample #v > 10000 #oneinth[k] #if the read count in this sample is > oneinten thr
                              lowabc += 1.0 if v > thrh[k] #add if the abundance is higher than 0.1% of the total sample count
                              #flat is the number of reads corresponding to 1% of the sample
                              flat = cos[k].to_f / 100.0 #* 5.0 #now a rep passes if its abundance is higher than 5% of the sample total
                              if v > cos[k].to_f / 100.0 #this is the 1% threshold in abundance in one sample
                                    fla1 = true  if v >= flat #keep if the read is 1% or more of reads in one sample
                                    countone += 1 #No of samples with this seq rep in >1% abundance
                              end
                              #override fla1 to 0.1% of the sample, if too few samples
                              #fla1 = true if v > cos[k].to_f / 10000.0 #keep all that are > 0.01%, useful for sets with only two samples
                              countzero +=1 if  v > cos[k].to_f / 1000.0 #No of samples with this seq rep in >0.1% abundance
                              zeroone = cos[k].to_f / 10000.0 #0.01% abundance threshold
                              countthr5 += 1 if v > zeroone #add to count thr5 
                        end
                  end
                  fla5 = true if countthr5 >= thr5 
                  lowfla = true if lowabc >= thr2 # at least n (thr2) samples with read abundance higher than 0.1% (thrh)
                  if (lowfla or fla5 or fla1)#or (tcc > mincom)) # n samples with abundance higher than thrh, OR present in more than thr5 samples OR read is 1% or more of all samples   
                        if (tcc > mincom and minreadnum)#mincom:thr on total copy No. tcc is total copy No. of this rep.
                              
                              keep.push(he)
                              repcps[he] = cps #seq rep -> hash with sample -> copy count
                              #removed header: #{he.gsub(/:|#/,"_")}
                              lin << "#{index},#{myd.faaha[he].seq},#{tcc.to_i},#{tsc.to_i},"
                              lcn << "#{index},#{myd.faaha[he].seq},#{tcc.to_i},#{tsc.to_i},"
                              faline << ">#{index}\n#{myd.faaha[he].seq}\n"   # collect the fasta in this line now
                              index += 1
                              cos.keys.sort.each do |m| #sample name
                                    if (cps[m].nil?)#||cps[m] < oneinth[m]) #if not there or below one in ten thousand copies
                                          lin << "0,"
                                          lcn << "0,"
                                    else
                                          lcn << "#{cps[m]}," #direct read count
                                          lin << "#{(cps[m]/cos[m])*100}," #percent of total reads in sample m 
                                    end
                              end
                              lin << "\n"
                              lcn << "\n"
                        end
                        
                  end
            end
            pretl = "Total read No,,,,,,,"
            cos.keys.sort.each{|k| pretl << "#{cos[k]},"} #add the reads per sample total to this line
            pretl << "\n"
            est = File.open("#{cle[0]}-table.csv","w")
            est.print lin
            est.close
            ess = File.open("#{cle[0]}-table-count.csv","w")
            ess.print lcn
            ess.close
            puts "retrieved #{lin.split(/\n/).length} reps"
            puts "retrieved #{index} seqs"
            lof.print("retrieved #{index} seqs\n")
            esf = File.open("#{cle[0]}-gutfilter.fa","w")
            esf.print faline
            esf.close
            lof.print "Min copy No: #{maxthr}\n"
            lof.close
      end
      
      def readspersample(uni) #the name file from the unique or pre.cluster output from mothur;  *.names
            stamp = Time.new.strftime("%Y%m%d-%H%M%S")
            #these are unique representatives and arrays of comma separated copies.
            coh = Hash.new(0.0) #counts the total number of seqs on each sample
            File.open(uni).each do |s|
                  re = s.chomp.split(/\s+/)[1]
                  re.split(/,/).each do |sa|
                        san = sa.split(/\||_/)[0]
                        coh[san] += 1.0 #counts the number of seqs per sample, gives also an array of non redundant sample names
                        #hal(reaa,he,san) #seq representative -> array of sample names
                  end
            end
            esc = File.open("failed-samples-#{stamp}.txt","w")
            esd = File.open("sample-read-Nos-#{stamp}.txt","w")
            esd.print "#{coh.keys.length} samples\nsample\treadNo\n"
            esc.print "sample\treadNo\n"
            coh.each do |k,v|
                  if v < 100 #low coverage test    
                        esc.print "#{k}\t#{v}\n"
                        coh.delete(k)
                  else
                        esd.print "#{k}\t#{v}\n"
                  end           
            end
            esc.close
            esd.close
            return coh
      end
      
      def unfiltered_abundance #reads the *names and gives abundance of these, for very low coverage samples
         #need array of all sample names.
         samplearray = `grep ">" all.fa | cut -f 1 -d "|" | sort |uniq `.gsub(/>/,"").split(/\n/)
         puts "X,#{samplearray.join(",")}\n" #header line
         saha = Hash.new
         samplearray.each{|k| saha[k] = 0} #initialize the hash with all samples as keys and 0 as default value
         File.open("all.trim.unique.precluster.names").each do |l|
            otu = l.split(/\s+/)[0] #the rep name
            l.split(/\s+/)[1].split(/,/).each do |sa|
                  san = sa.split(/\|/)[0]
                  saha[san] += 1
            end
            vals = "#{otu},"
            samplearray.each{|k| vals << "#{saha[k]},"}
            puts "#{vals}\n"
            samplearray.each{|k| saha[k] = 0} #re-initialize the hash to 0
         end
         #read the precluster name file
         #ditch the [0], split on /,/, for each 
      end
end
myg = Gut.new
myg.filter(ARGV[0],ARGV[1]) #preclustername file, precluster fasta file
#myg.readspersample(ARGV[0])
#myg.unfiltered_abundance
