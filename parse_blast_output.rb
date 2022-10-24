#!/usr/bin/ruby

require 'optparse'
require 'ostruct'


####################################################

## generate the HCV workflow as a shell script
full_command = ARGV.join(" ")
fastq_files = Array.new

#####################################################################
## parse arguments
class Parser
def self.parse(args)
	options = OpenStruct.new
	# default values
	options.type = "SILVA"
	options.consensus = "vote"
	options.weights = ""
	options.penaltyweight = 1.4
	options.useless_label_weight = 0.25
	options.full = false
	options.minscore = 0.002
	options.bestscore = 1e-199 # to prevent log(0)
	options.ignorestrain = false
	options.nlevels = 7
	options.minfrac = 0
	options.verbosity = 0

	# parse
	opt_parser = OptionParser.new do |opts|
		opts.banner = "Usage: parse_blast_output.rb [options] BLAST-FILE TAXONOMY-FILE"
		opts.on("-tTYPE", "--type=SILVA|UNITE|GTDB|NT", "[DEFAULT=SILVA] SILVA, UNITE, GTDB, or NCBI nt taxonomy file") do |val|
			options.type = val
		end
		opts.on("-cCONSENSUS", "--consensus=vote|weighted", "[DEFAULT=vote] Method to determine taxonomic assignment", "vote: classification with most votes wins", "weighted: weighted classification that includes suboptimal BLAST hits", "suboptimal hits are penalized by d = 1 - f*((Ex - Emin) / (Emax - Emin))", "where Ex is the log10(Evalue) of the current (suboptimal hit)", "Emin is the log10(Evalue) of the best hit", "Emax is the log10(value) of the worst hit", "and f is an optional constant (default=1)", "Weighted score is calculated as Wx*d", "where Wx is an optional taxa-specific weight (from --weight)") do |val|
			options.consensus = val
		end
		opts.on("-wWEIGHTS", "--weights", "File of weights to apply for taxonomic assignment (for plurality vote method)") do |val|
			options.weights = val
		end
		opts.on("-p", "--penaltyWeight PW", Float, "[DEFAULT=1.0] Constant modifier to suboptimal BLAST hit penalty function") do |val|
			options.penaltyweight = val
		end
		opts.on("--useless_label_weight ULW", Float, "[DEFAULT=0.25] Weight modifier for useless species-level taxonomic labels (e.g. Saccharomyces_sp)") do |val|
			options.useless_label_weight = val
		end
		opts.on("-f", "--full", "Print full results") do |val|
			options.full = val
		end
		opts.on("--ignore_strain", "Ignore strain when tabulating taxonomy") do |val|
			options.ignorestrain = val
			options.nlevels = 6
		end
		opts.on("--min_score MINSCORE", Float, "[DEFAULT=0.002] Minimum score to assign to suboptimal hits") do |val|
			options.minscore = val
		end
		opts.on("--best_score BESTSCORE", Float, "[DEFAULT=1e-199] Minimal (best scoring) evalue to accept (avoid zero)") do |val|
			options.bestscore = val
		end
		opts.on("--min_frac MINFRAC", Float, "[DEFAULT=0] Best hit must achieve at least this fraction of all votes/weights to be reported, otherwise report Ambiguous") do |val|
			options.minfrac = val
		end
		opts.on("-v", "--verbosity V", Integer, "[DEFAULT=0] 0=QUIET, 1=VERBOSE, 2=DEBUG") do |val|
			options.verbosity = val
		end
		opts.on("-h", "--help", "Prints this help") do
			puts opts
			exit
		end
	end
	opt_parser.parse!(args)
	return options
	end
end
options = Parser.parse(ARGV)
blast_fn, taxonomy_fn = ARGV

weights = Array.new
0.upto(options.nlevels) do |i|
	weights << Hash.new(1)
end

if options.weights != ""
	File.open(options.weights).each_line do |line|
		taxa, level, weight = line.chomp.split(/\t/)
		weights[level.to_i-1][taxa] = weight.to_f
#		puts "storing weight=#{weight.to_f} at weights[#{level.to_i-1}][#{taxa}]"
	end
end

####################################################

taxonomy = Hash.new
File.open(taxonomy_fn).each_line do |line|
	next unless line.chomp =~ /^>/
	if options.type == "SILVA" || options.type == "GTDB"
		arr = line.chomp.split(/\s+/)
		id = arr.shift.gsub(/^>/, "")
		tax = arr.join("_")
		taxonomy[id] = tax
	elsif options.type == "UNITE"
		#>UDB016649|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Thelephorales;f__Thelephoraceae;g__Thelephora;s__Thelephora_albomarginata|SH1502188.08FU
		arr = line.chomp.split(/\|/)
		id = arr[0].gsub(/^>/, "")
		tax = arr[1]
		tmp = tax.split(";").map { |lvlstr| lvlstr.split("__")[1]}
		tax = tmp.join(";")
		taxonomy[id] = tax
#		puts "storing taxonomy[#{id}] = #{tax}"
	elsif options.type == "NT"
		#M22247.1|B.fragilis ATP synthase beta-subunit gene, complete cds|817|Bacteroides fragilis
		arr = line.chomp.split(/\|/)
		id = arr[0]
		tax = arr[3]
		taxonomy[id] = tax
#		puts "storing taxonomy[#{id}] = #{tax}"
	end
end

## store BLAST hits
blast_best_scores = Hash.new(9999)
blast_worst_scores = Hash.new(-1)
blast_info = Hash.new # each element is an array holding the BLAST hits for that query
File.open(blast_fn).each_line do |line|
	#    'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', which is equivalent to the keyword 'std'
	query, subject, percent_id, len, num_mismatch, num_gapopen, qstart, qend, sstart, send, evalue, bitscore = line.chomp.split(/\t/)
	if evalue.to_f <= blast_best_scores[query]
		if evalue.to_f == 0
			evalue = options.bestscore
		end
		blast_best_scores[query] = evalue.to_f
	end
	if evalue.to_f > blast_worst_scores[query]
		blast_worst_scores[query] = evalue.to_f
	end
	if blast_info[query].nil?
		blast_info[query] = Array.new
	end
	blast_info[query] << line.chomp
end

## for each query sequence
blast_info.each_key { |qid|
	best_score = blast_best_scores[qid]
	worst_score = blast_worst_scores[qid]
#	puts "retrieving best_score=#{best_score} for qid=#{qid}"
	# voting scheme to determine consensus taxonomy at each level
	votes = Array.new
	0.upto(options.nlevels) do |lvl|
		votes << Hash.new(0) # a hash to store entries at each level
	end
	## for each BLAST hit
	blast_info[qid].each { |hit_str|
		# for 'vote' method, skip unless this is a best hit
		query, subject, percent_id, len, num_mismatch, num_gapopen, qstart, qend, sstart, send, evalue, bitscore = hit_str.chomp.split(/\t/)
		if evalue.to_f == 0
			evalue = options.bestscore
		end
		if (options.consensus=="vote")
			next unless evalue.to_f == best_score
		end
		if options.type == "SILVA" || options.type == "GTDB"
			taxstr = taxonomy[subject]
		else
			subject2 = subject.split(/\|/)[0]
			taxstr = taxonomy[subject2]
#			puts "subject=#{subject} subject2=#{subject2} taxstr=#{taxstr}"
		end
		i = 0
#		puts "taxstr=#{taxstr}"
		
		# compute weighted score if requested
		fd = 1
		if (options.consensus=="weighted" && evalue.to_f != best_score)
			fd = 1 - options.penaltyweight*(Math.log10(evalue.to_f) - Math.log10(best_score)) / (Math.log10(worst_score) - Math.log10(best_score))
			fd = fd
#			puts "fd=#{fd} evalue=#{evalue} options.penaltyweight=#{options.penaltyweight} best_score=#{best_score} worst_score=#{worst_score}"
			fd = [fd, options.minscore].max
			puts "calculated fd=#{fd} for query=#{query} taxstr=#{taxstr} evalue=#{evalue} best_score=#{best_score} worst_score=#{worst_score}" if options.verbosity==2
		end
		taxstr.split(";").each { |lvlstr|
#			puts "processing lvlstr=#{lvlstr} for taxstr=#{taxstr}"
			next if i > options.nlevels # happens for Eukaryota strings (e.g. Eukaryota;Archaeplastida;Chloroplastida;Charophyta;Phragmoplastophyta;Streptophyta;Embryophyta;Tracheophyta;Spermatophyta;Magnoliophyta;Solanales;Solanum;Solanum_lycopersicum_(tomato))
			next if lvlstr=="uncultured_bacterium" || lvlstr=="uncultured_organism" || lvlstr=="microorganism" || lvlstr=="metagenome" || lvlstr=="uncultured" || lvlstr=="unidentified" || lvlstr=="unclassified" # skip useless classification hits
			modifier_useless_label = 1 # optional modifier for useless taxonomic labels (e.g. Saccharomyces_sp)
			if (i==7 && options.ignorestrain)
				puts "ignoring strain for lvlstr=#{lvlstr}" if options.verbosity==2
				lvlstr = lvlstr.split("_")[0..1].join("_")
				puts "new lvlstr=#{lvlstr}" if options.verbosity==2
			end
			if (lvlstr.split("_")[1] == "sp" || lvlstr.split("_")[1] == "sp.")
				modifier_useless_label = options.useless_label_weight
			end
			votes[i][lvlstr] += 1*(weights[i][lvlstr])*fd*modifier_useless_label
			puts "adding vote #{lvlstr} at level #{i} with value #{1*(weights[i][lvlstr])*fd*modifier_useless_label}" if options.verbosity==2
			i+=1
		}
	}
	
	## get best consensus taxonomy for this query sequence
	tax_consensus = Array.new
	0.upto(options.nlevels) do |lvl|
		consensus = Array.new
		maxweight = 0
		totalweight = 0
		votes[lvl].each { |k, v|
			totalweight = totalweight + v
			if ((!options.full) && v == votes[lvl].values.max)
				consensus << k
				maxweight = maxweight + v
#				puts "found max value #{v} at key=#{k} value=#{v} level=#{lvl}"
				puts "votes[lvl]=#{votes[lvl].values} with lvl=#{votes[lvl].keys}" if options.verbosity==2
			elsif (options.full)
				consensus << "#{k}[#{v}]"
			end
		}
		# check against options.minfrac
		if !options.full
#			puts "trying to divide maxweight=#{maxweight} totalweight=#{totalweight} for lvl=#{lvl}"
			if maxweight > 0 && totalweight > 0 && (1.0*maxweight / totalweight) < options.minfrac
#				puts "got maxweight=#{maxweight} totalweight=#{totalweight} frac=#{1.0*maxweight/totalweight} for lvl=#{lvl}"
				consensus.clear
				consensus << "Ambiguous"
			end
		end
		tax_consensus << consensus.join("+")
	end
	taxstr = tax_consensus.join(";")
	eval = best_score
	puts "#{qid}\t#{taxstr}\t#{eval}"
}
