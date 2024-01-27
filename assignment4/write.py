with open('potential2.txt', 'w') as file:
    # Write a hundred zeros to the file, each on a new line
    x_R = 1
    x_L = -1.00000000000000000000000
    N = 100
    dx = (x_R - x_L) / (N - 1)
    
    for _ in range(100):
        # q = (-1/(x_L)) * 1000
        p = 1/2*(x_L**2)*1000
        x_L += dx
        # use q for hydrogen atom (change x_L to be slightly above zero, but not zero), or use p for harmonic oscillator
        # file.write(f"{p}\n")
    
    # to write a 100 zeros
    for _ in range(100):
        file.write("0.0\n")
    
    
    # to write a tunnel thingy 
    # for _ in range(40):
    #     file.write("0.0\n")
    
    # for _ in range(20):
    #     file.write("10.0\n")
        
    # for _ in range(40):
    #     file.write("0.0\n")
        
    
