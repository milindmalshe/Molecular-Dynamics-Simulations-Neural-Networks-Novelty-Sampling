
global Terminate
Terminate = 0;
hProg = waitbar(0, 'Simulating', 'CreateCancelBtn', 'Stop' );
waitbar( 0.0 );
while 0 == Terminate

for i=1:10000
    for j=1:10000
        for k=1:10000
            a(i,j,k)=1;
        end
    end
    i
end

end
delete( hProg );
Terminate = 0;
