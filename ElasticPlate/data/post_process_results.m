root = 'N=500/'

isx = importdata(strcat(root,'imported_sx.txt'));
rsx = importdata(strcat(root,'result_sx.txt'));
isy = importdata(strcat(root,'imported_sy.txt'));
rsy = importdata(strcat(root,'result_sy.txt'));
isxy = importdata(strcat(root,'imported_sxy.txt'));
rsxy = importdata(strcat(root,'result_sxy.txt'));

dx = isx(:,5) - rsx(:,5);
dy = isy(:,5) - rsy(:,5);
dxy = isxy(:,5) - rsxy(:,5);
rcss = [isx(:,2)';isx(:,3)';dx';dy';dxy']';
csvwrite(strcat(root,'ss_diff.csv'),rcss);