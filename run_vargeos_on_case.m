

 init=load_untouch_nii('init.nii');
 t1=load_untouch_nii('t1.nii');
 t1ce=load_untouch_nii('t1ce.nii');
 

 init.img=~(init.img>0); 

 Gamma=double((t1ce.img-t1.img)>1);
 Gamma=(~Gamma);
    
 rex=vargeos(t1ce.img,init.img,Gamma);


 t1ce.img=rex;
 save_untouch_nii(t1ce,'vargeos.nii'])
