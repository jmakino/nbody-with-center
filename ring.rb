#
# place m rings of n particles with given mass in witdth dr
include Math
m=20
n=1000
dr=0.1
mass=5e-8
if ARGV.length == 4 then
  m=ARGV[0].to_i
  n=ARGV[1].to_i
  dr=ARGV[2].to_f
  mass=ARGV[3].to_f
end
STDERR.print "m=#{m} n=#{n} dr=#{dr} mass=#{mass}\n"
def generate_ring(r, n, mass, idoffset)
  dtheta=Math::PI*2.0/n
  vcoef = 1.0/sqrt(r)
  n.times{|i|
    theta = dtheta*i;
    print idoffset+i," ",mass, " ", r*cos(theta), " ", r*sin(theta), " ",0, " "
    print  vcoef*sin(theta), " ", -vcoef*cos(theta)," ", 0, "\n"
  }
end

dx = dr/(m-1)
print "0\n"
print n*m, "\n"
m.times{|i|
  r = 1.0+dx*i;
  id=i*n
  generate_ring(r, n, mass, id)
}

    
  

