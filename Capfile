###
# Basic Limma analysis on the data, Annotation, Ontology analysis etc.

require 'catpaws'


# Your AWS credentials
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']

# ec2 url to use - this is set in ~/.bash_profile to the eu-west. 
set :ec2_url, ENV['EC2_URL']

# The type of instances you wish to start
set :ami, 'ami-7e5c690a'  #EC2 eu-west-1 32bit Maverick
set :instance_type, 'm1.small'

# The group you want to start these instances in. Under most circumstances this
# should be a unique name that doesn't clash with group names of already running 
# instances.
set :group_name, 'ns5_vs_ns5dastro_lumixpn'

# How many instances do you want?
# Note - most of the EBS: tasks won't work with more that one instance yet. I'm working on it ;)
set :nhosts, 1

# Where to start your instances.  Note that this needs to be a zone in your region 
# and if you change ec2_url you'll need to change this accordingly. 
# You don't *have* to define it, but if you're using an EBS vol you need to make
# sure it's in the same zone as your instances, and the easiest way to do that is
# to specify the zone you're using here.
set :availability_zone, 'eu-west-1a' 

# Connection options for your instances
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']   # key *name* 
set :key_file, ENV['EC2_KEYFILE']  # path to keyfile


# If you don't want to use your home directory, specify a working dir here.
# If you aren't using an EBS-backed image, your instance storage will probably 
# be accessible under /mnt and you'll run out of room if you start dumping 
# files in your home dir. 
set :working_dir, '/mnt/work'



# If defined, EBS:create will create a volume from this snapshot
set :snap_id, `cat SNAPID`.chomp 

# If defined, EBS:attach will attach this volume
set :vol_id, `cat VOLUMEID`.chomp

# How big (in GB) would you like your EBS volume (smaller snapshots will be resized)
set :ebs_size, 1 

# Where would you like your EBS volume attached?
set :dev, '/dev/sdf'

# And where would you like it mounted?
set :mount_point, '/mnt/data'



#make a new EBS volume from this snap 
#cap EBS:create

#and mount your EBS
#cap EBS:attach
#cap EBS:format_xfs
#cap EBS:mount_xfs


desc "a test task"
task :a_test, :roles => group_name do
  results = capture("ls /home")
  puts results
end
before 'a_test', 'EC2:start'


#if you want to keep the results

#cap EBS:snapshot

#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




