import imageio

new_dir = "iteration_2000"


maxIter = 1000 
plots = 50

frames = []
t = maxIter//plots

print("\nMaking Gif ...")
for i in range(t):
    num = "{0:0=3d}".format(i)
    print(num)
#     image = imageio.v2.imread(f"./" + new_dir + "/fluid." + num + ".png")
#     frames.append(image)

# imageio.mimsave("./" + new_dir + ".gif", frames,duration = 80)