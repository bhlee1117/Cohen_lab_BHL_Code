function [regist_pat] = register_pattern_test(L,W,offsets)
% L = 1024; W = 768;
t_r = 10;% registration pattern thickness
t_m = 7; % marking number thickness
l_r = 50;% registration pattern length
l_m = 20;% marking number length
w_r = 50;% registration pattern width
w_m = 10;% marking number width
off_r_x = offsets(1); %registration pattern x offset;
off_r_y = offsets(2); %registration pattern y offset;
off_m_x = 30; %registration pattern x offset;
off_m_y = 40; %registration pattern y offset;
sep_m=25;
one_tip=5;
off_one = 6;


pat = false(L,W);
roi00_sq1_x = [1:50]+off_r_x; roi00_sq0_x = [10:50]+off_r_x; 
roi00_sq1_y = [1:50]+off_r_y; roi00_sq0_y = [10:50]+off_r_y; 
roi00_mark1_x = [find((abs((1:W)-off_m_x-off_r_x))<10) find((abs((1:W)-(off_m_x+sep_m)-off_r_x))<10)]; 
roi00_mark0_x = [find((abs((1:W)-off_m_x-off_r_x))<3) find((abs((1:W)-(off_m_x+sep_m)-off_r_x))<3)];
roi00_mark1_y = [find((abs((1:L)-off_m_y-off_r_y))<20) find((abs((1:L)-off_m_y-off_r_y))<20)]; 
roi00_mark0_y = [find((abs((1:L)-off_m_y-off_r_y))<13) find((abs((1:L)-off_m_y-off_r_y))<13)];

roi01_sq1_x = [1:50]+off_r_x; roi01_sq0_x = [10:50]+off_r_x; 
roi01_sq1_y = L-[0:49]-off_r_y; roi01_sq0_y = L-[10:49]-off_r_y; 
roi01_mark1_x0 = find((abs((1:W)-off_m_x-off_r_x))<10);      roi01_mark1_x1 = [find((abs((1:W)-(off_m_x+sep_m)-off_r_x))<8)]; 
roi01_mark0_x0 = find((abs((1:W)-off_m_x-off_r_x))<3);       roi01_mark0_x1= find((abs((1:W)-(off_m_x+sep_m-one_tip)-off_r_x))<3);
roi01_mark1_y0 = find((abs((1:L)-(L-off_m_y-off_r_y)))<20);      roi01_mark1_y1 = [find((abs((1:L)-(L-off_m_y-off_r_y)))<20)]; 
roi01_mark0_y0 = find((abs((1:L)-(L-off_m_y-off_r_y)))<13);      roi01_mark0_y1 = find((abs((1:L)-(L-(off_m_y-off_one)-off_r_y)))<14);

roi10_sq1_x = W-[0:49]-off_r_x; roi10_sq0_x = W-[10:49]-off_r_x; 
roi10_sq1_y = [1:50]+off_r_y; roi10_sq0_y = [10:50]+off_r_y; 
roi10_mark1_x0 = find((abs((1:W)-(W-off_m_x-off_r_x)))<10);      roi10_mark1_x1 = find((abs((1:W)-(W-(off_m_x+sep_m)-off_r_x)))<8); 
roi10_mark0_x0 = find((abs((1:W)-(W-off_m_x-off_r_x)))<3);       roi10_mark0_x1= find((abs((1:W)-(W-(off_m_x+sep_m+one_tip)-off_r_x)))<3);
roi10_mark1_y0 = find((abs((1:L)-off_m_y-off_r_y))<20);      roi10_mark1_y1 = find((abs((1:L)-off_m_y-off_r_y))<20); 
roi10_mark0_y0 = find((abs((1:L)-off_m_y-off_r_y))<13);      roi10_mark0_y1 = find((abs((1:L)-(off_m_y+off_one)-off_r_y))<14);

pat00 = pat; pat00_mark = pat;
pat01 = pat; pat01_mark = pat; 
pat10 = pat; pat10_mark = pat; 

pat00(roi00_sq1_y,roi00_sq1_x)=true; 
pat00(roi00_sq0_y,roi00_sq0_x)=false; 
pat00_mark(roi00_mark1_y,roi00_mark1_x)=true; 
pat00_mark(roi00_mark0_y,roi00_mark0_x)=false;

pat01(roi01_sq1_y,roi01_sq1_x)=true; 
pat01(roi01_sq0_y,roi01_sq0_x)=false; 
pat01_mark(roi01_mark1_y0,roi01_mark1_x0)=true; 
pat01_mark(roi01_mark0_y0,roi01_mark0_x0)=false;
pat01_mark(roi01_mark1_y1,roi01_mark1_x1)=true; 
pat01_mark(roi01_mark0_y1,roi01_mark0_x1)=false;


pat10(roi10_sq1_y,roi10_sq1_x)=true; 
pat10(roi10_sq0_y,roi10_sq0_x)=false; 
pat10_mark(roi10_mark1_y0,roi10_mark1_x0)=true; 
pat10_mark(roi10_mark0_y0,roi10_mark0_x0)=false;
pat10_mark(roi10_mark1_y1,roi10_mark1_x1)=true; 
pat10_mark(roi10_mark0_y1,roi10_mark0_x1)=false;

regist_pat = pat00|pat00_mark|pat01|pat01_mark|pat10|pat10_mark;
% figure
% imagesc(pat00|pat00_mark|pat01|pat01_mark|pat10|pat10_mark)
% axis equal

end