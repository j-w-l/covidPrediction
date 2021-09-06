% Jonathan Lee, PS5.

% *PURPOSE: Constructs 2x2 contact matrix. This should (theoretically) be
% proportional to the WAIFW matrix of beta coefficients. (WAIFW_M = Pr[infection] * Contact_M).
% *DEPENDENCIES: us_contact.csv from S1 Dataset in Prem et al. '17: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697.
% *CAVEAT: Prem et al.'s original contact matrix was 8x8. Our 2x2
% generalization may be reducitve, or even wrong.

% Logistics.
filename='us_contact.csv';
opts=detectImportOptions(filename);
contact_matrix=readtable(filename, opts);

% Computing entry C11 of contact matrix.
entries = 0;
C11 = 0;
for i=1:13
    avg = 0;
    for j=2:14
        avg = avg + contact_matrix{i, j};
        entries = entries + 1;
    end
    C11 = C11 + (avg./13);
end
C11 = C11/13;

% Computing entry C12.
entries = 0;
C12 = 0;
for i=1:13
    avg = 0;
    for j=15:17
        avg = avg + contact_matrix{i, j};
        entries = entries + 1;
    end
    C12 = C12 + (avg./3);
end
C12 = C12/13;

% Computing entry C21.
entries = 0;
C21 = 0;
for i=14:16
    avg = 0;
    for j=2:14
        avg = avg + contact_matrix{i, j};
        entries = entries + 1;
    end
    C21 = C21 + (avg./13);
end
C21 = C21/3;

% Computing entry C22.
entries = 0;
C22 = 0;
for i=14:16
    avg = 0;
    for j=15:17
        avg = avg + contact_matrix{i, j};
        entries = entries + 1;
    end
    C22 = C22 + (avg./3);
end
C22 = C22/3;

% To print the computed matrix entries:
[C11, C12, C21, C22]
% writetable(contact_matrix, 'adjusted_contact_matrix.csv')