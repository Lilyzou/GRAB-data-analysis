%This file will take data from an individual session and repackage it into
%a separate data structure so it can be compared with another session from
%a different animal. Make sure to load the collected load file first.

%%Please copy this first
session_count=0;

%% Copy and paste this section when you want to load in a session
session_count = session_count+1;
Simon_Session_Struct
disp('Please input the session number (order # in the load file; e.g. in 0122-0123-0124.mat if you want 0123, please type 2)');
session_number = input(sprintf('Session number %d:', session_count));
extracted_session(session_count) = session(session_number);

%huh okay that was less than I thought I needed