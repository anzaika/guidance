#! /usr/bin/ruby
#/usr/bin/env ruby

tree = gets.strip.gsub( /_.*?:/, ":" ).gsub(/:0,/, ":0.00,").gsub(/:0\)/, ":0.00\)").gsub(/[0-9]\.[0-9]*e-[0-9]*/, "0.0").gsub(/[0-9]*e-[0-9]*/, "0.0")
# Buggy line:
#tree = gets.strip.gsub( /_.*?:/, ":" ).gsub(/:0,/, ":0.00,").gsub(/:0\)/, ":0.00\)").gsub(/[0-9]\.[0-9]*e-[0-9][0-9]/, "0.0")


#puts "tree = " +  tree

memi = [-1,-1]
leni = [-1,-1]

c = 0
while tree.index( /\(/ ) 

	tree.sub!( /\(([0-9]*?):(\-?[0-9]*?\.[0-9]*?),([0-9]*?):(\-?[0-9]*?\.[0-9]*?)\)/, "XXX" )
	memi[0] = $1.to_i
	leni[0] = $2.to_f
	memi[1] = $3.to_i
	leni[1] = $4.to_f

	if leni[0] > 3 || leni[1] > 3 then
		STDERR.puts "Please check the scale of branch length!"
		STDERR.puts "1PAM must be 0.01, not 1"
#		exit 1
	end

#	puts "subtree = " + $&

	if memi[1] < memi[0] then
		memi.reverse!
		leni.reverse!
	end

	tree.sub!( /XXX/, memi[0].to_s )
#	puts "tree = " + tree

	c = c+1
	printf( "%5d %5d %10.5f %10.5f\n", memi[0], memi[1], leni[0], leni[1] )
	if c > 1000000 then
	   print( "Endless loop bug!!!  printed > 1000000 lines\n" )
	   exit 1
	end

end

