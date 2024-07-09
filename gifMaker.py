import imageio


directory = "Flow_no_roll_flags_obstacle_6000_it"
maxIter = 6000
plots = 50

frames = []
t = maxIter//plots

print("Making Gif ...")
for i in range(t):
    # print(i)
    num = "{0:0=5d}".format(i)
    image = imageio.v2.imread(f"./" + directory + "/fluid_" + str(num) + ".png")
    frames.append(image)

imageio.mimsave("./" + directory + "/_" + directory + ".gif", frames, duration = 80)

print("\nDone.")