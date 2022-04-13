% function to find the indices of propagator elements. This function is
% used when the rodrigues() function is used to calulate propagator and
% propagator derivative elements for use in ESCALADE and QOALA methods.
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu

function indices=propagator_ind(space,basis,n_spins)

indices=cell(1,n_spins);

N_row=kron([1 1 1],[1 2 3]);
N_col=kron([1 2 3],[1 1 1]);
I_dim=kron([1 1 1],[1 1 1]);

I=[1 1 1 1];
N=[0 1 2 3];

if strcmp(space,'liouville') && strcmp(basis,'sphten')

    E=eye(n_spins);
    for n=1:n_spins % spins
        prop_ind=1;
        N_mat=E([n_spins:-1:n_spins-n+2 1:n_spins-n+1],:);
        for m1=1:size(N_mat,1) % sums
            kron_rows=1;
            kron_cols=1;
            for m2=1:size(N_mat,2) % krons
                if N_mat(m1,m2)
                    if m2==1
                        row_cell=N_row;
                        col_cell=N_col;
                    else
                        row_cell=N;
                        col_cell=N;
                    end
                else
                    if m2==1
                        row_cell=I_dim;
                        col_cell=I_dim;
                    else
                        row_cell=I;
                        col_cell=I;
                    end
                end
                kron_rows=kron(row_cell,kron_rows);
                kron_cols=kron(col_cell,kron_cols);
            end

            prop_ind=prop_ind+[kron_rows(:) kron_cols(:)].*(4^(n_spins-m1));
        end
        indices{n}=prop_ind;
    end

    if n_spins==1
        indices{1}=[1 1; indices{1}];
    end

elseif strcmp(space,'liouville') && strcmp(basis,'zeeman')

    row_index= kron(I,[1 2 3 4]);
    col_index= kron([1 2 3 4],I);

    prop_ind=[row_index(:) col_index(:)];

    if n_spins==1
        indices{1}=prop_ind;
    elseif n_spins==2

        row_index= kron(I,kron(I,N))*4^1+...
            kron([1 2 3 4]     ,kron(I,I))*4^0;

        col_index= kron(I,kron(N,I))*4^1+...
            kron([1 2 3 4]     ,kron(I,I))*4^0;

        indices{1}=[row_index(:) col_index(:)];

        row_index= kron(I,kron(I,[1 2 3 4]))*4^0+...
            kron(N,kron(I,I))*4^1;

        col_index= kron(I,kron([1 2 3 4],I))*4^0+...
            kron(N,kron(I,I))*4^1;

        indices{2}=[row_index(:) col_index(:)];

    elseif n_spins==3


        row_index= kron(I,kron(I,N_row))*4^2+...
            kron(I,kron(N,I_dim))*4^0+...
            kron(N,kron(I,I_dim))*4^1;

        col_index= kron(I,kron(I,N_col))*4^2+...
            kron(I,kron(N,I_dim))*4^0+...
            kron(N,kron(I,I_dim))*4^1;

        indices{1}=[row_index(:) col_index(:)];

        row_index= kron(I,kron(I,N_row))*4^1+...
            kron(I,kron(N,I_dim))*4^0+...
            kron(N,kron(I,I_dim))*4^2;

        col_index= kron(I,kron(I,N_col))*4^1+...
            kron(I,kron(N,I_dim))*4^0+...
            kron(N,kron(I,I_dim))*4^2;

        indices{2}=[row_index(:) col_index(:)];

        row_index= kron(I,kron(I,N_row))*4^0+...
            kron(I,kron(N,I_dim))*4^1+...
            kron(N,kron(I,I_dim))*4^2;

        col_index= kron(I,kron(I,N_col))*4^0+...
            kron(I,kron(N,I_dim))*4^1+...
            kron(N,kron(I,I_dim))*4^2;

        indices{3}=[row_index(:) col_index(:)];

    elseif n_spins==4


        row_index= kron(I,kron(I,kron(I,N_row)))*4^3+...
            kron(I,kron(I,kron(N,I_dim)))*4^1+...
            kron(I,kron(N,kron(I,I_dim)))*4^2+...
            kron(N,kron(I,kron(I,I_dim)))*4^0;

        col_index= kron(I,kron(I,kron(I,N_col)))*4^3+...
            kron(I,kron(I,kron(N,I_dim)))*4^1+...
            kron(I,kron(N,kron(I,I_dim)))*4^2+...
            kron(N,kron(I,kron(I,I_dim)))*4^0;

        indices{1}=[row_index(:) col_index(:)];

        row_index= kron(I,kron(I,kron(I,N_row)))*4^2+...
            kron(I,kron(I,kron(N,I_dim)))*4^1+...
            kron(I,kron(N,kron(I,I_dim)))*4^3+...
            kron(N,kron(I,kron(I,I_dim)))*4^0;

        col_index= kron(I,kron(I,kron(I,N_col)))*4^2+...
            kron(I,kron(I,kron(N,I_dim)))*4^1+...
            kron(I,kron(N,kron(I,I_dim)))*4^3+...
            kron(N,kron(I,kron(I,I_dim)))*4^0;

        indices{2}=[row_index(:) col_index(:)];

        row_index= kron(I,kron(I,kron(I,N_row)))*4^1+...
            kron(I,kron(I,kron(N,I_dim)))*4^0+...
            kron(I,kron(N,kron(I,I_dim)))*4^2+...
            kron(N,kron(I,kron(I,I_dim)))*4^3;

        col_index= kron(I,kron(I,kron(I,N_col)))*4^1+...
            kron(I,kron(I,kron(N,I_dim)))*4^0+...
            kron(I,kron(N,kron(I,I_dim)))*4^2+...
            kron(N,kron(I,kron(I,I_dim)))*4^3;

        indices{3}=[row_index(:) col_index(:)];

        row_index= kron(I,kron(I,kron(I,N_row)))*4^0+...
            kron(I,kron(I,kron(N,I_dim)))*4^1+...
            kron(I,kron(N,kron(I,I_dim)))*4^2+...
            kron(N,kron(I,kron(I,I_dim)))*4^3;

        col_index= kron(I,kron(I,kron(I,N_col)))*4^0+...
            kron(I,kron(I,kron(N,I_dim)))*4^1+...
            kron(I,kron(N,kron(I,I_dim)))*4^2+...
            kron(N,kron(I,kron(I,I_dim)))*4^3;

        indices{4}=[row_index(:) col_index(:)];

    else
        error('zeeman-liouv not yet coded for n>1 spins')
    end

elseif strcmp(space,'hilbert') && strcmp(basis,'zeeman')

    row_index= kron([1 1],[1 2]);
    col_index= kron([1 2],[1 1]);

    prop_ind=[row_index(:) col_index(:)];

    if n_spins==1
        indices{1}=prop_ind;
    elseif n_spins==2

        row_index= kron([1 1]     ,kron([1 1],[0 1]))*2^1+...
            kron([1 2]     ,kron([1 1],[1 1]))*2^0;

        col_index= kron([1 2]     ,kron([1 1],[1 1]))*2^0+...
            kron([1 1]     ,kron([0 1],[1 1]))*2^1;
        %
        %    (1,1)     0.000000000000000 + 0.001899221753683i y (1,1)
        %    (3,1)     0.000022741896039 - 0.499995189852507i y (3,1)
        %    (2,2)    -0.000022741896039 - 0.499995189852507i n (1,3)
        %    (4,2)     0.000000000000000 - 0.001899221753683i n (3,3)
        %    (1,3)     0.000000000000000 + 0.001899221753683i n (2,2)
        %    (3,3)     0.000022741896039 - 0.499995189852507i n (4,2)
        %    (2,4)    -0.000022741896039 - 0.499995189852507i y (2,4)
        %    (4,4)     0.000000000000000 - 0.001899221753683i y (4,4)

        indices{1}=[row_index(:) col_index(:)];

        row_index= kron([1 1]     ,kron([1 1],[1 2]))*2^0+...
            kron([0 1]     ,kron([1 1],[1 1]))*2^1;

        col_index= kron([1 1]     ,kron([1 2],[1 1]))*2^0+...
            kron([0 1]     ,kron([1 1],[1 1]))*2^1;

        indices{2}=[row_index(:) col_index(:)];

    else
        %warning('zeeman-hilb not yet coded for n>2 spins')
    end

else

    error(['propagator_ind not coded for the ' basis ' basis in a ' space ' space formalism'])

end
end
