

% define platten 

jojo =Platten('A','N',[1,0,0],[0,0,1])


%jojo2 =Platten('A','N',[1,0,0],[0,0,1],0.2,0.1)


myTrans=Transducers([1,2,3,4,5],['s','s','s','r','r'],[0,1,2,3,4],['A','A','A','A','A'],[1,2,3,4,5]);

Block=struct;
Block.L_E=25;
Block.L_N=25;
Block.L_T=25;

xyz = myTrans.calculate_global_coordinates({jojo},Block)


