Toro225

if(x > 0.3) {
    u[0] = 0.125;
    u[1] = 0;
    u[2] = 0.1 / 0.4;
} else {
    u[0] = 1.;
    u[1] = 0.75;
    u[2] = 1./ 0.4+0.5*0.75*0.75;
}