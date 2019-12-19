#
# place m rings of n particles with given mass in witdth dr
include Math
m=20
n=1000
dr=0.1
mass=5e-8
vrandscale=0.0
if ARGV.length >= 4 then
  m=ARGV[0].to_i
  n=ARGV[1].to_i
  dr=ARGV[2].to_f
  mass=ARGV[3].to_f
  vrandscale=ARGV[4].to_f if ARGV.length >= 5

end
STDERR.print "m=#{m} n=#{n} dr=#{dr} mass=#{mass} vrandscale=#{vrandscale}\n"
def scaled_rand(scale)
  Random::rand(1.0)*scale*2 -scale
end
def generate_ring(r, n, mass, idoffset, vrandscale)
  dtheta=Math::PI*2.0/n
  vcoef = 1.0/sqrt(r)
  n.times{|i|
    theta = dtheta*i;
    print idoffset+i," ",mass, " ", r*cos(theta), " ", r*sin(theta), " ",0, " "
    print  vcoef*sin(theta)+scaled_rand(vrandscale), " ",
           -vcoef*cos(theta)+scaled_rand(vrandscale)," ", scaled_rand(vrandscale), "\n"
  }
end

dx = dr/(m-1)
Random::srand(12345)
print "0\n"
print n*m, "\n"
m.times{|i|
  r = 1.0+dx*i;
  id=i*n
  generate_ring(r, n, mass, id, vrandscale)
}

    
  

