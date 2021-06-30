imshow(im1);
hold on;
plot(m_pt_1(1:8,1),m_pt_1(1:8,2),'go')

homo_x2 = [m_pt_2(1:8,:), ones(8,1,'double')];

epiLines2 = [ homo_x2*HF ];
points = lineToBorderPoints(epiLines2,size(im1));

%Plot the epipolar lines on the image
line(points(:,[1,3])',points(:,[2,4])');
truesize;

