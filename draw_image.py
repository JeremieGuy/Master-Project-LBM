
from matplotlib import image 
from matplotlib import pyplot as plt 
  
# to read the image stored in the working directory 
data = image.imread('./Monitoring/u_roll_301_151_new_collision_2000/system.png') 

line1_1 = [150,150]
line1_2 = [1,52]

plt.plot(line1_1,line1_2,color="red", linewidth=3)
plt.imshow(data)
plt.show()