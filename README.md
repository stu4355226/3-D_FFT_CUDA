# Multithreading-Analysis-for-3-D-Fast-Fourier-Transform
A 3D-FFT implementation on C and CUDA.
Author    
-------------
CHIA-CHE, LEE   
Computer Engineering, Florida institute of technology.   
Email: stu4355226@gmail.com    
    
Target   
-------------
Platform: Visual Studio 2013.
Language: C, CUDA.  
   
Program Flow   
-------------
   
![program_flow](/images/gram_flow.jpg)   
   
This project is using CUDA and C languages to implement 3D Fast Fourier Transform and 3D inverse Fast Fourier Transform.     
The goal is compare the executuon speed and anylysis the advantages and disadvantages on using CUDA for multthreading.    
In 3D-FFT folder, the ouput file is "fft3rm.d" file.     
In the first exeution, the program will generate a 3D array with ramdom numbers.    
         
![original](/images/original.jpg)      
     
Then add the transfer data right after the original data.     
     
![transfer](/images/transfer.jpg)      
      
      
Then add the inverse data  after the transfer data.    
   
![inverse](/images/inverse.jpg)   
   
if the incerse data is the same as the original data, we can say the program is executing correctly.   
   
to check the data for CUDA 3D-FFT.   
you need to open file"fft3_cuda.d", the file is made in the same way of C 3D-FFT does.   

Execution 
-------------
   
![execution](/images/execution.jpg)   
   
Result   
-------------
   
![Result](/images/result.jpg) 
   
The result is showing that in the small matrix computation, CUDA is slower becase of moving data from DDR to GDDR taking too much time.    
It's on my expection. If we don't count the time that moves data from DDR to GDDR, CUDA is faster in any cases.    
