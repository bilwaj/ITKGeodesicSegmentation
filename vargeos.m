function[result]=vargeos(Inp,mask,Gamma)
%geos(Inp,mask,Gamma)
Gamma=double(Gamma);
Inp=double(Inp);
if sum(mask(:)==0)+sum(mask(:)==1)~=max(size(mask(:)))
 sprintf('Mask does not contain just 0s and 1s. Please check the mask and try again')
 result='Error'
else
  result=10000000*mask;     
  [Z,Y,X]=size(result)
  C_f_arr=zeros(14,1);C_b_arr=zeros(14,1);
  for z=2:Z-1
      for y=2:Y-1
          for x=2:X-1
                                        C_f_arr(14)=result(z,y,x);
					C_f_arr(1)=result(z-1,y-0,x-0)+sqrt(1.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y-0,x-0))*(Inp(z,y,x)-Inp(z-1,y-0,x-0))));
					C_f_arr(2)=result(z-0,y-1,x-0)+sqrt(1.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y-1,x-0))*(Inp(z,y,x)-Inp(z-0,y-1,x-0))));
					C_f_arr(3)=result(z-0,y-0,x-1)+sqrt(1.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y-0,x-1))*(Inp(z,y,x)-Inp(z-0,y-0,x-1))));
					C_f_arr(4)=result(z-1,y-1,x-0)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y-1,x-0))*(Inp(z,y,x)-Inp(z-1,y-1,x-0))));
					C_f_arr(5)=result(z-1,y-0,x-1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y-0,x-1))*(Inp(z,y,x)-Inp(z-1,y-0,x-1))));
					C_f_arr(6)=result(z-0,y-1,x-1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y-1,x-1))*(Inp(z,y,x)-Inp(z-0,y-1,x-1))));
					C_f_arr(7)=result(z-1,y-1,x-1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y-1,x-1))*(Inp(z,y,x)-Inp(z-1,y-1,x-1))));
					C_f_arr(8)=result(z-1,y+1,x-0)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y+1,x-0))*(Inp(z,y,x)-Inp(z-1,y+1,x-0))));
					C_f_arr(9)=result(z-1,y+1,x-1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y+1,x-1))*(Inp(z,y,x)-Inp(z-1,y+1,x-1))));
					C_f_arr(10)=result(z-0,y+1,x-1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y+1,x-1))*(Inp(z,y,x)-Inp(z-0,y+1,x-1))));
					C_f_arr(11)=result(z+1,y+1,x-1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y+1,x-1))*(Inp(z,y,x)-Inp(z+1,y+1,x-1))));
					C_f_arr(12)=result(z+1,y-0,x-1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y-0,x-1))*(Inp(z,y,x)-Inp(z+1,y-0,x-1))));
					C_f_arr(13)=result(z+1,y-1,x-1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y-1,x-1))*(Inp(z,y,x)-Inp(z+1,y-1,x-1))));
					result(z,y,x)=min(C_f_arr);
		
          end
      end
  end
  for z=Z-2:-1:2
      for y=Y-2:-1:2
          for x=X-2:-1:2
                 
				        C_b_arr(14)=result(z,y,x);
					C_b_arr(1)=result(z+1,y-0,x-0)+sqrt(1.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y-0,x-0))*(Inp(z,y,x)-Inp(z+1,y-0,x-0))));
					C_b_arr(2)=result(z-0,y+1,x-0)+sqrt(1.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y+1,x-0))*(Inp(z,y,x)-Inp(z-0,y+1,x-0))));
					C_b_arr(3)=result(z-0,y-0,x+1)+sqrt(1.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y-0,x+1))*(Inp(z,y,x)-Inp(z-0,y-0,x+1))));
					C_b_arr(4)=result(z+1,y+1,x-0)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y+1,x-0))*(Inp(z,y,x)-Inp(z+1,y+1,x-0))));
					C_b_arr(5)=result(z+1,y-0,x+1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y-0,x+1))*(Inp(z,y,x)-Inp(z+1,y-0,x+1))));
					C_b_arr(6)=result(z-0,y+1,x+1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y+1,x+1))*(Inp(z,y,x)-Inp(z-0,y+1,x+1))));
					C_b_arr(7)=result(z+1,y+1,x+1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y+1,x+1))*(Inp(z,y,x)-Inp(z+1,y+1,x+1))));
					C_b_arr(8)=result(z+1,y-1,x-0)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y-1,x-0))*(Inp(z,y,x)-Inp(z+1,y-1,x-0))));
					C_b_arr(9)=result(z+1,y-1,x+1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z+1,y-1,x+1))*(Inp(z,y,x)-Inp(z+1,y-1,x+1))));
					C_b_arr(10)=result(z-0,y-1,x+1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-0,y-1,x+1))*(Inp(z,y,x)-Inp(z-0,y-1,x+1))));
					C_b_arr(11)=result(z-1,y-1,x+1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y-1,x+1))*(Inp(z,y,x)-Inp(z-1,y-1,x+1))));
					C_b_arr(12)=result(z-1,y-0,x+1)+sqrt(2.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y-0,x+1))*(Inp(z,y,x)-Inp(z-1,y-0,x+1))));
					C_b_arr(13)=result(z-1,y+1,x+1)+sqrt(3.0+Gamma(z,y,x)*((Inp(z,y,x)-Inp(z-1,y+1,x+1))*(Inp(z,y,x)-Inp(z-1,y+1,x+1))));
					result(z,y,x)=min(C_b_arr);
          end
      end
  end
  MAX_VAL=result(1,1,1);
  result(result==MAX_VAL)=0;
  
end
