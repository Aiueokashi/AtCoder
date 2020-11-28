
DD ax,ay,bx,by,cx,cy;
in(ax,ay,bx,by,cx,cy);
out(Point(ax,ay)-Point(bx,by));
out(simple_ccw(Point(ax,ay),Point(bx,by),Point(cx,cy)));
out(Apporonius(Point(ax,ay),Point(bx,by),1,2));

//-------------------------------------------------------//
DD ax,ay,bx,by,cx,cy;
in(ax,ay,bx,by,cx,cy);
vector<Point> vec(2);
rep(i,2){
cin >> vec[i].x >> vec[i].y ;
}
out(Point(ax,ay)-Point(bx,by));
out(simple_ccw(Point(ax,ay),Point(bx,by),Point(cx,cy)));
out(ccw(Point(ax,ay),Point(bx,by),Point(cx,cy)));
out(Apporonius(Point(ax,ay),Point(bx,by),1,2));
cout << Line(vec);
