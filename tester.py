import Camera
C = Camera.Camera(testpattern=True, subarray=100)
I = C.ccds[0]
B = C.cartographer
for i in B.possibleinputs:
    B.point(0,0,i)
    for o in B.possibleoutputs:
        B.quote(o)
I.expose(write=True)
