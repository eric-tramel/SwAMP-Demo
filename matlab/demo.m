addpath bin common examples

fprintf('Welcome to SwAMP demo script!\n')
setup;

% Standard projector
fprintf(['First, we''ll try both SwAMP and AMP with Gaussian i.i.d. ' ...
    'projectors of zero mean and 1/N variance.\n'])
fprintf('(Press any key.)\n')
pause;

fprintf('\n')
fprintf('o Executing demo for an i.i.d. zero mean projector...\n')
nzm(0);
fprintf('\n')

% Non-zero mean projector
fprintf(['As we can see, both work and reach the same MSE.\n'])
fprintf(['Let''s now see what happens if we set the mean to 20 / N...\n'])
fprintf('(Press any key.)\n')
pause;

fprintf('\n')
fprintf('o Executing demo for an i.i.d. *positive mean* projector...\n')
nzm(20);
fprintf('\n')

% Low-rank projector
fprintf('AMP diverges in this case, while SwAMP remains converging to a good MSE.\n')
fprintf(['What about low-rank sensing matrices? Let''s use a projector F = P * Q, ' ...
    'following the paper''s prescription and with \\eta < \\alpha.\n'])
fprintf('(Press any key.)\n')
pause;

fprintf('\n')
fprintf('o Executing demo for a highly correlated, low-rank projector...\n')
lr(0.5);
fprintf('\n')

% 1-bit CS
fprintf('Once again, SwAMP shows to be effective, while AMP doesn''t work.\n')
fprintf(['We now move to 1-bit CS, starting with a consistency check: GAMP ' ...
    'and G-SwAMP should return similar results for Gaussian i.i.d. matrices ' ...
    'of zero mean.\n'])
fprintf('(Press any key.)\n')
pause;

fprintf('\n')
fprintf('o Executing 1-bit CS demo for an i.i.d. zero mean projector...\n')
bit(0);
fprintf('\n')

fprintf('But what about positive mean matrices?\n')
fprintf('(Press any key.)\n')
pause;

fprintf('\n')
fprintf('o Executing 1-bit CS demo for an i.i.d. *positive mean* projector...\n')
bit(20);
fprintf('\n')

% Sparse matrices
fprintf(['Finally, we show that SwAMP can be almost as fast as AMP for sparse matrices. ' ...
    'Let''s consider a Gaussian i.i.d. matrix with 20%% of non-zero elements; ' ...
    'look at the runtimes!\n'])
fprintf('(Press any key.)\n')
pause;

fprintf('\n')
fprintf('o Executing demo for sparse projector with 20%% of non-zero elements...\n')
sprs(0.2)
fprintf('\n')

fprintf('Thanks for your attention! :-)\n')
