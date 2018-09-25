function prediction(infile_name)
num_identifiers = 1;
infile = fopen(['../../' infile_name], 'r');  
while ~feof(infile) % continue to iterate until file pointer 'infile' is at the end of the file
	file_line = fgetl(infile);
	if strcmp(file_line(1), '>')
		space_loc = strfind(file_line, ' ');
		if ~isempty(space_loc)
  			allidentifiers{num_identifiers} = file_line(2:space_loc(1)-1);
  		else
  			allidentifiers{num_identifiers} = file_line(2:end);
		end
  		num_identifiers = num_identifiers+1;
	end
end
fclose(infile);

%%%%%%%
%%%%%%%
%%%%long-range
%%%%%%%
%%%%%%%

if ~exist('../results/long')
	mkdir('../results/long');
end

if ~exist('../results/long_40')
	mkdir('../results/long_40');
end

load('./Final_model/long_model_40.mat');
load('./Final_model/max_min_vals_long_40_50.mat');

for i = 1:num_identifiers-1
	identifier = allidentifiers{i};
	
	resultoutfile = fopen(['../results/long_40/' identifier '.DeepRCon'], 'w');

	Dtest = load(['../features/long/' identifier '_features']);
	
	if ~isempty(Dtest)
		 
		for c = 91:290
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val1)/(max_val1-min_val1);
		end
		for c = 413:472
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val2(c-412))/(max_val2(c-412)-min_val2(c-412));
		end

		infile_obj = fopen(['../features/long/' identifier '_features_obj']);
		C = textscan(infile_obj, '%s %s'); 
		fclose(infile_obj);
		fprintf(resultoutfile, '#res1 res2      contact_prob      no-contact_prob\n');
		N = size(Dtest, 1); 

		Dtest = [Dtest ones(N,1)];
		w1probs = 1./(1 + exp(-Dtest*w1)); w1probs = [w1probs  ones(N,1)];
		w2probs = 1./(1 + exp(-w1probs*w2));
		clear w1probs;
		w2probs = [w2probs ones(N,1)];
		w3probs = 1./(1 + exp(-w2probs*w3)); 
		clear w2probs;
		w3probs = [w3probs  ones(N,1)];
		targetsout = exp(w3probs*w_class);
		clear w3probs;
		targetsout = targetsout./repmat(sum(targetsout,2),1,2);


		for j = 1:size(targetsout, 1)  
			fprintf(resultoutfile, '%s %s      %f      %f\n', C{1}{j}, C{2}{j}, targetsout(j, 1), targetsout(j, 2));
		end
	end
		
	fclose(resultoutfile);
end
clear w1 w2 w3 w_class max_val2 max_val1 min_val1 min_val2;

if ~exist('../results/long_50')
	mkdir('../results/long_50');
end

load('./Final_model/long_model_50.mat');
load('./Final_model/max_min_vals_long_40_50.mat');

for i = 1:num_identifiers-1
	identifier = allidentifiers{i};
	
	resultoutfile = fopen(['../results/long_50/' identifier '.DeepRCon'], 'w');

	Dtest = load(['../features/long/' identifier '_features']);
	
	if ~isempty(Dtest)
		 
		for c = 91:290
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val1)/(max_val1-min_val1);
		end
		for c = 413:472
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val2(c-412))/(max_val2(c-412)-min_val2(c-412));
		end

		infile_obj = fopen(['../features/long/' identifier '_features_obj']);
		C = textscan(infile_obj, '%s %s'); 
		fclose(infile_obj);
		fprintf(resultoutfile, '#res1 res2      contact_prob      no-contact_prob\n');
		N = size(Dtest, 1); 

		Dtest = [Dtest ones(N,1)];
		w1probs = 1./(1 + exp(-Dtest*w1)); w1probs = [w1probs  ones(N,1)];
		w2probs = 1./(1 + exp(-w1probs*w2));
		clear w1probs;
		w2probs = [w2probs ones(N,1)];
		w3probs = 1./(1 + exp(-w2probs*w3)); 
		clear w2probs;
		w3probs = [w3probs  ones(N,1)];
		targetsout = exp(w3probs*w_class);
		clear w3probs;
		targetsout = targetsout./repmat(sum(targetsout,2),1,2);


		for j = 1:size(targetsout, 1)  
			fprintf(resultoutfile, '%s %s      %f      %f\n', C{1}{j}, C{2}{j}, targetsout(j, 1), targetsout(j, 2));
		end
	end
		
	fclose(resultoutfile);
end
clear w1 w2 w3 w_class max_val2 max_val1 min_val1 min_val2;


if ~exist('../results/long_60')
	mkdir('../results/long_60');
end

load('./Final_model/long_model_60.mat');
load('./Final_model/max_min_vals_long_60.mat');

for i = 1:num_identifiers-1
	identifier = allidentifiers{i};
	
	resultoutfile = fopen(['../results/long_60/' identifier '.DeepRCon'], 'w');

	Dtest = load(['../features/long/' identifier '_features']);
	
	if ~isempty(Dtest)
		 
		for c = 91:290
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val1)/(max_val1-min_val1);
		end
		for c = 413:472
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val2(c-412))/(max_val2(c-412)-min_val2(c-412));
		end

		infile_obj = fopen(['../features/long/' identifier '_features_obj']);
		C = textscan(infile_obj, '%s %s'); 
		fclose(infile_obj);
		fprintf(resultoutfile, '#res1 res2      contact_prob      no-contact_prob\n');
		N = size(Dtest, 1); 

		Dtest = [Dtest ones(N,1)];
		w1probs = 1./(1 + exp(-Dtest*w1)); w1probs = [w1probs  ones(N,1)];
		w2probs = 1./(1 + exp(-w1probs*w2));
		clear w1probs;
		w2probs = [w2probs ones(N,1)];
		w3probs = 1./(1 + exp(-w2probs*w3)); 
		clear w2probs;
		w3probs = [w3probs  ones(N,1)];
		targetsout = exp(w3probs*w_class);
		clear w3probs;
		targetsout = targetsout./repmat(sum(targetsout,2),1,2);


		for j = 1:size(targetsout, 1)  
			fprintf(resultoutfile, '%s %s      %f      %f\n', C{1}{j}, C{2}{j}, targetsout(j, 1), targetsout(j, 2));
		end
	end
		
	fclose(resultoutfile);
end
clear w1 w2 w3 w_class max_val2 max_val1 min_val1 min_val2;


%%%%%%%
%%%%%%%
%%%%medium-range
%%%%%%%
%%%%%%%

if ~exist('../results/medium')
	mkdir('../results/medium');
end

load('./Final_model/medium_model.mat');
load('./Final_model/max_min_vals_medium.mat');

for i = 1:num_identifiers-1
	identifier = allidentifiers{i};
	
	resultoutfile = fopen(['../results/medium/' identifier '.DeepRCon'], 'w');

	Dtest = load(['../features/medium/' identifier '_features']);
	
	if ~isempty(Dtest)
		 
		for c = 91:290
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val1)/(max_val1-min_val1);
		end
		for c = 415:474
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val2(c-414))/(max_val2(c-414)-min_val2(c-414));
		end

		infile_obj = fopen(['../features/medium/' identifier '_features_obj']);
		C = textscan(infile_obj, '%s %s'); 
		fclose(infile_obj);
		fprintf(resultoutfile, '#res1 res2      contact_prob      no-contact_prob\n');
		N = size(Dtest, 1); 

		Dtest = [Dtest ones(N,1)];
		w1probs = 1./(1 + exp(-Dtest*w1)); w1probs = [w1probs  ones(N,1)];
		w2probs = 1./(1 + exp(-w1probs*w2)); w2probs = [w2probs ones(N,1)];
		w3probs = 1./(1 + exp(-w2probs*w3)); w3probs = [w3probs  ones(N,1)];
		targetsout = exp(w3probs*w_class);   
		targetsout = targetsout./repmat(sum(targetsout,2),1,2);


		for j = 1:size(targetsout, 1)  
			fprintf(resultoutfile, '%s %s      %f      %f\n', C{1}{j}, C{2}{j}, targetsout(j, 1), targetsout(j, 2));
		end
	end
		
	fclose(resultoutfile);
end
clear w1 w2 w3 w_class max_val2 max_val1 min_val1 min_val2;

%%%%%%%
%%%%%%%
%%%%short-range
%%%%%%%
%%%%%%%

if ~exist('../results/short')
	mkdir('../results/short');
end

load('./Final_model/short_model.mat');
load('./Final_model/max_min_vals_short.mat');

for i = 1:num_identifiers-1
	identifier = allidentifiers{i};
	
	resultoutfile = fopen(['../results/short/' identifier '.DeepRCon'], 'w');

	Dtest = load(['../features/short/' identifier '_features']);
	
	if ~isempty(Dtest)
		 
		for c = 51:250
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val1)/(max_val1-min_val1);
		end
		for c = 351:410
			Dtest(:,c) = 1.0*(Dtest(:,c)-min_val2(c-350))/(max_val2(c-350)-min_val2(c-350));
		end

		infile_obj = fopen(['../features/short/' identifier '_features_obj']);
		C = textscan(infile_obj, '%s %s'); 
		fclose(infile_obj);
		fprintf(resultoutfile, '#res1 res2      contact_prob      no-contact_prob\n');
		N = size(Dtest, 1); 

		Dtest = [Dtest ones(N,1)];
		w1probs = 1./(1 + exp(-Dtest*w1)); w1probs = [w1probs  ones(N,1)];
		w2probs = 1./(1 + exp(-w1probs*w2)); w2probs = [w2probs ones(N,1)];
		w3probs = 1./(1 + exp(-w2probs*w3)); w3probs = [w3probs  ones(N,1)];
		targetsout = exp(w3probs*w_class); 
		targetsout = targetsout./repmat(sum(targetsout,2),1,2);


		for j = 1:size(targetsout, 1)  
			fprintf(resultoutfile, '%s %s      %f      %f\n', C{1}{j}, C{2}{j}, targetsout(j, 1), targetsout(j, 2));
		end
	end
		
	fclose(resultoutfile);
end


quit;
