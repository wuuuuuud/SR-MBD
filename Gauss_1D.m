function ker = Gauss_1D( HWHH )
%���ݸ����İ�߿��һ�룬����һ����˹ģ���ˣ�����һ��
%input:  HWHH - half width of the half high
%output: ker - the result of gaussian bulry kernel

xx = [-round(3*HWHH):1:round(3*HWHH)];%��ɢ��������֤�ǶԳƵ�
ker = 1 / sqrt(2 * pi * HWHH^2) * exp(( -xx.^2)/(2 * HWHH^2 ));
nor = sum(sum( ker ));
ker = ker / nor;
end
