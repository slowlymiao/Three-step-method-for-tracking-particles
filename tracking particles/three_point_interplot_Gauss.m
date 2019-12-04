function [xm]=three_point_interplot_Gauss(col,r_1,r,r1,x_win)
    xm=col+(log(r_1)-log(r1))/(log(r_1)+log(r1)-2*log(r))/2-x_win;
 
end