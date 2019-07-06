function bst_events = utest_get_test_bst_events()
% Return a BST event struct with corner cases:
%   - multiple conditions
%   - some have empty notes / channels tags
%   - some have non-empty notes with weird chars
%   - some selected, some not
%
% Maintain this function everytime bst event structure changes
% -> to add new fields use expected default values so that backward compatibility 
%    is guarenteed. test_load_save_backward_compatibility should yell at
%    you if not.
bst_events = db_template('event');
bst_events(1).label = 'condition1';
bst_events(1).color = [0., 0., 0.];
bst_events(1).epochs = [1];
bst_events(1).times = [10.3;10.9];
bst_events(1).reactTimes = [];
bst_events(1).select = 1;
bst_events(1).channels = {{}};
bst_events(1).notes = {''};

bst_events(2).label = 'condition2';
bst_events(2).color = [0., 1., 0.3];
bst_events(2).epochs = [1 1 1];
bst_events(2).times = [1.3 7.5 12.6];
bst_events(2).reactTimes = [];
bst_events(2).select = 0;
bst_events(2).channels = {{'S1D1', 'S2D3'}, {}, {'S121', 'S4D3'}};
bst_events(2).notes = {sprintf('\t!:/\\,#%%*+{}[]wazaéàè?`''&$'), ...
                   '', '12345678910-=;,1,2||'};

end