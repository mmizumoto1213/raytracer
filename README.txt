(1) Matthew Mizumoto. I am working alone
(2) Username: mmizumoto
    Email: mmizumoto@ucsd.edu
(3) https://raviucsdgroup.s3.amazonaws.com/hw3/6bc10199c7cbd8964e42ed1694b75e1c/20240320015223/index.html
(4) The code was ran using Visual Studio Code. It was compiled using "g++ main.cpp readfile.cpp Transform.cpp -o main.exe -L. -lFreeImage -O3" In my readfile I never included a case to read the output so every time I executed my code I manually named all of my scenes using the third argument in argv. (ex of execution to name a file "./main.exe scene7.test scene7.png")
(5)
(6) I worked alone was not able to implement any acceleration structure. However, I did attempt to create one (bounding boxes) which failed and is commented out but still left in my code. Without the acceleration structure I was able to run my scene7 in roughly 40 minutes.