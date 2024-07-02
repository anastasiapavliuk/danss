zs = []
ys = []
xs = []
for i in range(1,3):
    f = open("./sumevents/coord_result/coord_2/_"+str(i)+".txt", 'r')
    c = 0
    for line in f:
        c+=1
        if (c > 3 and line[1] != "*" and line[0] == "*"):
            spl_line = line.split("*")
            zs.append(spl_line[-2])
            ys.append(spl_line[-3])
            xs.append(spl_line[-4])

f = open("final_koord.txt", 'w')
for i in range(len(xs)):
   f.write(xs[i] + " " + ys[i] + " " + zs[i] + '\n')
