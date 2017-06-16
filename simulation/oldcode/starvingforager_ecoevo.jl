function starvingforager_ecoevo(L,dim,initsize,t_term,alpha,K,sigma_mean,sigma_sd,rho_mean,rho_sd,m,lambda,mu)
  #Read in packages/function
  #ipbc :: torus movement
  #include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/ipbc.jl")


  #Initiate values
  #Lattice dimension

  S = (L-2)^dim;

  #intisize after rounding
  r_initsize = Int(round(initsize/3));

  #Initial state values
  # Find_vec = zeros(Int64,r_initsize) + 2;
  # Hind_vec = zeros(Int64,r_initsize) + 1;
  # Rind_vec = zeros(Int64,r_initsize);

  time_out = Array(Float64,0);
  N_out = Array(Int64,0);

  NF = r_initsize;
  NH = r_initsize;
  NR = r_initsize;

  N = NF + NH + NR;

  #Initial densities
  F = NF/S;
  H = NH/S;
  R = NR/S;
  prop_out = Array{Float64}(3,1);
  prop_out[:,1] = [F,H,R];
  push!(N_out,N);



  #Initial location values
  #2/3 of initsize will be foragers; 1/3 will be resources
  #replace = false because want only one resource per site

  Floc_vec = sample(collect(1:S),r_initsize);
  Hloc_vec = sample(collect(1:S),r_initsize);
  #Use replace=false to ensure that no location is chosen twice for resources only
  Rloc_vec = sample(collect(1:S),r_initsize,replace=false);

  #Define mutation variance (this might depend on the rate values!)
  mutation = 0.001;

  #DEFINE TRAIT VEC
  #Keep track of F and H trait_vec... maybe use later.
  if (sigma_sd > 0)
    d_sigma = Normal(sigma_mean,sigma_sd);
    F_sigma_vec = rand(d_sigma,r_initsize);
    H_sigma_vec = rand(d_sigma,r_initsize);
    d_mutation = Normal(0,mutation);
  else
    F_sigma_vec = ones(r_initsize)*sigma_mean;
    H_sigma_vec = ones(r_initsize)*sigma_mean;
    d_mutation = Bernoulli(0);
  end
  if (rho_sd > 0)
    d_rho = Normal(rho_mean,rho_sd);
    F_rho_vec = rand(d_rho,r_initsize);
    H_rho_vec = rand(d_rho,r_initsize);
    d_mutation = Normal(0,mutation);
  else
    F_rho_vec = ones(r_initsize)*rho_mean;
    H_sigma_vec = ones(r_initsize)*rho_mean;
    d_mutation = Bernoulli(0);
  end

  #Set initial sigma and rho
  sigma = mean([F_sigma_vec;H_sigma_vec]);
  rho = mean([F_rho_vec;H_rho_vec]);

  trait_out = Array{Float64}(2,1);
  trait_out[:,1] = [sigma,rho];


  push!(time_out, 0);

  #NEED TO ENSURE THAT RESOURCES ARE PLACED 1 PER SITE
  #HAVE A noresourcesites vector that accounts for sites WITHOUT resources
  #Updated each timestep
  noresourcesites = collect(1:S);
  #Deletes locations of sites with resources from the list of open sites (noresourcesites)
  deleteat!(noresourcesites,sort(Rloc_vec));

  t = 0;
  next_time_int = t + 1;


  #push!(prop_out,copy(prop));

  tic = 0;

  # F_pr_line = 0.0;
  # H_pr_line = 0.0;
  # R_pr_line = 0.0;
  # id = Array(Float64,1);

  while t < (t_term-1)
    tic = tic + 1;

    # #ERRORS
    # if length(ind_vec) < 1
    #   print("Extinction has occured at t=",round(t,2)," and loop ",tic)
    #   break
    # end

    #Calculate Rate
    Rate = F*(lambda + sigma*(K-R)) + H*(rho*R + mu) + R*(alpha*(K-R) + (rho*H + m*F)); # (1 - (N/S)) +

    # TESTING
    # Rate = (NF/N)*(lambda + sigma*(K-R)) + (NH/N)*(rho*R + mu) + (NR/N)*(alpha*(K-R) + (F+H));

    dt = 1/(Rate*N);
    if Rate == 0
      println("Welcome to Daisy World")
      break
    end



    #Randomly select an individual (R,S,F) with probability 1/N
    #ind thus represents the POSITION of the individual
    #Update total number of individuals
    pr_line = zeros(6);
    #Update the total
    #Events
    pr_line[1] = (lambda*F)/Rate;
    #2  Starvation
    pr_line[2] = pr_line[1] + (sigma*(K-R)*F)/Rate;
    #3  Recruitment
    pr_line[3] = pr_line[2] + (rho*H*R)/Rate;
    #4  Death
    pr_line[4] = pr_line[3] + (mu*H)/Rate;
    #5  Resource Growth
    pr_line[5] = pr_line[4] + (alpha*R*(K-R))/Rate;
    #6  Resource consumption
    pr_line[6] = pr_line[5] + (m*F*R)/Rate;

    draw_event = rand();

    #1  Consumer Reproduction (F)
    if draw_event < pr_line[1]
      #Add individual to end of Find_vec
      # push!(Find_vec,2);
      #Add random location to Floc_vec
      location = rand(collect(1:S));
      push!(Floc_vec,location);

      #Introduce new mutated traits
      id = rand(collect(1:NF)); #This is the individual that is reproducing
      new_sigma = copy(F_sigma_vec[id]) + rand(d_mutation);
      push!(F_sigma_vec,new_sigma);
      new_rho = copy(F_rho_vec[id]) + rand(d_mutation);
      push!(F_rho_vec,new_rho);

      #Update
      NF = NF + 1;
    end

    #2  Starvation (F)
    if draw_event >= pr_line[1] && draw_event < pr_line[2]
      #Randomly draw F position from Find_vec
      id = rand(collect(1:NF));
      location = Floc_vec[id];
      sigma_id = copy(F_sigma_vec[id]);
      rho_id = copy(F_rho_vec[id]);
      #Delete that individual from Find_vec/Floc_vec; add to Hind_vec/Hloc_vec
      # deleteat!(Find_vec,id);
      deleteat!(Floc_vec,id);
      # push!(Hind_vec,1);
      push!(Hloc_vec,location);

      #Track trait values
      deleteat!(F_sigma_vec,id);
      push!(H_sigma_vec,sigma_id);
      deleteat!(F_rho_vec,id);
      push!(H_rho_vec,rho_id);

      #Update
      NF = NF - 1;
      NH = NH + 1;
    end

    #3  Recruitment/Recover (H)
    if draw_event >= pr_line[2] && draw_event < pr_line[3]
      #Randomly draw H position from Hind_vec
      id = rand(collect(1:NH));
      location = Hloc_vec[id];
      sigma_id = copy(H_sigma_vec[id]);
      rho_id = copy(H_rho_vec[id]);
      #Delete that individual from Hind_vec/Hloc_vec; add to Find_vec/Floc_vec
      # deleteat!(Hind_vec,id);
      deleteat!(Hloc_vec,id);
      # push!(Find_vec,2);
      push!(Floc_vec,location);

      #Track trait values
      deleteat!(H_sigma_vec,id);
      push!(F_sigma_vec,sigma_id);
      deleteat!(H_rho_vec,id);
      push!(F_rho_vec,rho_id);

      #Remove a resource unit
      idR = rand(collect(1:NR));
      locationR = Rloc_vec[idR];
      #Delete from vector
      # deleteat!(Rind_vec,id);
      deleteat!(Rloc_vec,idR);
      push!(noresourcesites,locationR);

      #Update
      NH = NH - 1;
      NF = NF + 1;
      NR = NR - 1;
    end

    #4  Death (H)
    if draw_event >= pr_line[3] && draw_event < pr_line[4]
      #Randomly draw H position
      id = rand(collect(1:NH));
      #Delete from vector
      # deleteat!(Hind_vec,id);
      deleteat!(Hloc_vec,id);
      deleteat!(H_sigma_vec,id);
      deleteat!(H_rho_vec,id);

      #Update
      NH = NH - 1;
    end

    #5  Resource Growth (R)
    if draw_event >= pr_line[4] && draw_event < pr_line[5]
      #Draw a random location WITHOUT A RESOURCE
      newresourcepos = rand(collect(1:length(noresourcesites)));
      location = noresourcesites[newresourcepos];

      # push!(Rind_vec,0);
      push!(Rloc_vec,location);

      #Update
      NR = NR + 1;
      #Delete the filled resource postion from list of noresourcesites
      deleteat!(noresourcesites,newresourcepos);
    end

    #6  Resource consumption (R)
    if draw_event >= pr_line[5] && draw_event < pr_line[6]

      #Draw a random R postion
      id = rand(collect(1:NR));
      location = Rloc_vec[id];
      #Delete from vector
      # deleteat!(Rind_vec,id);
      deleteat!(Rloc_vec,id);

      #Update
      NR = NR - 1;
      #Update noresourcesites (add site that resource was eliminated from)
      push!(noresourcesites,location);

    end


    #Recalculate the total # of the foragers and resources
    N = NF + NH + NR;

    #Recalculate the densities of each
    F = NF/S;
    H = NH/S;
    R = NR/S;

    prop = [F,H,R];

    #push!(prop_out,copy(prop));
    prop_out = hcat(prop_out,prop);

    #Recalculate the average sigma and rho
    sigma = mean([F_sigma_vec;H_sigma_vec]);
    rho = mean([F_rho_vec;H_rho_vec]);

    trait = [sigma,rho];
    trait_out = hcat(trait_out,trait);

    #Advance time
    t = t + dt;

    if t > next_time_int
      println("time= ",round(t,2))
      next_time_int = round(t+1,0)
    end

    #Make sure only values are pushed to the output
    # ind_vec_new = copy(ind_vec);
    # loc_vec_new = copy(loc_vec);

    #Update output
    # push!(ind_out,ind_vec_new);
    # push!(loc_out,loc_vec_new);
    push!(time_out,t);
    push!(N_out,N);

    #ERRORS
    #Break loop if extinction occurs
    if N < 1
      println("Extinction has occured at t=",round(t,2)," and loop ",tic)
      break
    end
    #Break loop if resources go extinct and the other populations run away
    if NR == 0 && (NF + NH) > S*2
      println("Runaway growth has occured at t=",round(t,2)," and loop ",tic)
      break
    end



  end #end while loop over t

  println("Simulation successful at t= ",round(t,2)," and loop= ",tic)

  # return ind_out,loc_out,time_out,prop_out,N_out
  return time_out,prop_out,N_out,trait_out
end #end function
