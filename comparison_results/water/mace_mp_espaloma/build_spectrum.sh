BROAD=0.0025


grep Inten $1 | cut -c 16-150 | gawk '{printf "%12.4f\n%12.4f\n%12.4f\n", $1, $2, $3}' > intensity.dat

grep Frequ $1 | cut -c 16-150 | gawk '{printf "%12.4f\n%12.4f\n%12.4f\n", $1, $2, $3}' > frequencies.dat



paste frequencies.dat intensity.dat  > lines.dat


gawk '{

        count = 0;
        for (i=10; i < 4000; i = i + 0.25)
        { 

           argument=  -1.0 * '$BROAD' *  (i - $1)^2;

            #printf "debug %d    %lf    %lf\n", i, argument, 0.05*(i - $1)^2;	   
	   
	   if (argument > -745) 
           {

             I[count] = I[count] + $2 * exp(  - '$BROAD' * ( i - $1)^2)
	     
	     #printf "debug %d    %lf    %lf\n", i, I[count], exp(  -1.0 * '$BROAD' * ( i - $1)^2);
	   }  

           count++;  
        }
	
      }; END {

        count = 0;
        for (i=10; i < 4000; i = i + 0.25)
        { 
          printf ("%8.4lf    %16.12e\n", i, I[count]); 
	  count++;  
        }
      }' lines.dat
  
