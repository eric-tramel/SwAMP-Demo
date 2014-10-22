function setup
    newline = sprintf('\n');
    src = regexprep(ls('src/*/*.c'), '[ \t\n]+', ' ');
    cxx = 'mex';
    cflags = '-g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99"';

    fprintf('Compiling SwAMP... ')
    
    eval([cxx ' ' cflags ' -output bin/swamp src/swamp.c ' src]);
    eval([cxx ' ' cflags ' -output bin/swgamp src/swgamp.c ' src]);   
    
    % We should probably check and see if we have any lingering output
    % files, just in case. This seems to be an issue with building OSX.
    eval('delete *.o');
    
fprintf('Done!\n')
end
