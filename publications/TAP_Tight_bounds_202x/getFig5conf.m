function conf = getFig5conf(ID)

switch ID
    case 1
        conf.M = 8;
        conf.w = taylorwin(conf.M,3,-20); conf.w = conf.w / sum(conf.w);
        conf.theta = deg2rad(40);
        conf.ampErr = 0.1;
        conf.phaErr = deg2rad(5);
        conf.mtlCpl = 0.01;
        conf.xL = [-0.3 0.3];
        conf.yL = [-0.155 0.155];
end

end