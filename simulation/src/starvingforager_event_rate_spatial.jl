function starvingforager_event_rate_spatial(L,dim,initsize,t_term,alpha,K,sigma,rho,m,lambda,mu,Df,Dh,write_out)
  #Read in packages/function
  #ipbc :: torus movement
  include("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/ipbc.jl")


  #Initiate values
  #Lattice dimension

  S = (L+2)^dim;
  S_func = L^dim;
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
  F = NF/S_func;
  H = NH/S_func;
  R = NR/S_func;
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

  #Record the location of F/H/R in this vector, where each element is a site.
  #F and H can have more than one element in a site; R can only have 1.
  Fsite_vec = zeros(S);
  Hsite_vec = zeros(S);
  for i=1:r_initsize
    Fsite_vec[Floc_vec[i]] += 1
    Hsite_vec[Hloc_vec[i]] += 1
  end

  Rsite_vec = zeros(S);
  Rsite_vec[Rloc_vec] = 1;

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
    Rate = F*(lambda + sigma*(K-R) + Df) + H*(rho*R + mu + Dh) + R*(alpha*(K-R) + (rho*H + m*F)); # (1 - (N/S)) +

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
    pr_line = zeros(8);
    #Update the total
    #Events
    #1 F growth
    pr_line[1] = (lambda*F)/Rate;
    #2  Starvation
    pr_line[2] = pr_line[1] + (sigma*(K-R)*F)/Rate;
    #3  F Movement
    pr_line[3] = pr_line[2] + (Df*F)/Rate;
    #3  Recruitment
    pr_line[4] = pr_line[3] + (rho*H*R)/Rate;
    #4  Death
    pr_line[5] = pr_line[4] + (mu*H)/Rate;
    #5  H Movement
    pr_line[6] = pr_line[5] + (Dh*H)/Rate;
    #5  Resource Growth
    pr_line[7] = pr_line[6] + (alpha*R*(K-R))/Rate;
    #6  Resource consumption
    pr_line[8] = pr_line[7] + ((rho*H + m*F)*R)/Rate;

    draw_event = rand();

    #1  Consumer Reproduction (F)
    if draw_event < pr_line[1]
      #Add individual to end of Find_vec
      # push!(Find_vec,2);
      #Add random location to Floc_vec
      location = rand(collect(1:S));
      push!(Floc_vec,location);

      #update site vectors
      Fsite_vec[location] += 1;

      #Update
      NF = NF + 1;
    end

    #2  Starvation (F)
    if draw_event >= pr_line[1] && draw_event < pr_line[2]
      #Randomly draw F position from Find_vec
      id = rand(collect(1:NF));
      location = Floc_vec[id];
      #Delete that individual from Find_vec/Floc_vec; add to Hind_vec/Hloc_vec
      # deleteat!(Find_vec,id);
      deleteat!(Floc_vec,id);
      # push!(Hind_vec,1);
      push!(Hloc_vec,location);

      #update site vectors
      Fsite_vec[location] -= 1;
      Hsite_vec[location] += 1;

      #Update
      NF = NF - 1;
      NH = NH + 1;
    end

    #3 F Movement
    if draw_event >= pr_line[2] && draw_event < pr_line[3]
      #Randomly draw F position
      id = rand(collect(1:NF));
      location = Floc_vec[id];

      nn = ipbc(location,L);
      new_location = rand(nn);
      Floc_vec[id] = new_location;

      #update site vectors
      Fsite_vec[location] -= 1;
      Fsite_vec[new_location] += 1;
    end

    #4  Recruitment/Recover (H)
    if draw_event >= pr_line[3] && draw_event < pr_line[4]
      #Randomly draw H position from Hind_vec
      id = rand(collect(1:NH));
      location = Hloc_vec[id];
      #Delete that individual from Hind_vec/Hloc_vec; add to Find_vec/Floc_vec
      # deleteat!(Hind_vec,id);
      deleteat!(Hloc_vec,id);
      # push!(Find_vec,2);
      push!(Floc_vec,location);

      Fsite_vec[location] += 1;
      Hsite_vec[location] -= 1;

      #Update
      NH = NH - 1;
      NF = NF + 1;
    end

    #5  Death (H)
    if draw_event >= pr_line[4] && draw_event < pr_line[5]
      #Randomly draw H position
      id = rand(collect(1:NH));
      location = Hloc_vec[id];
      #Delete from vector
      # deleteat!(Hind_vec,id);
      deleteat!(Hloc_vec,id);

      #update site vectors
      Hsite_vec[location] -= 1;

      #Update
      NH = NH - 1;
    end

    #6  H Movement
    if draw_event >= pr_line[5] && draw_event < pr_line[6]
      #Randomly draw H position
      id = rand(collect(1:NH));
      location = Hloc_vec[id];

      nn = ipbc(location,L);
      new_location = rand(nn);
      Hloc_vec[id] = new_location;

      #update site vectors
      Hsite_vec[location] -= 1;
      Hsite_vec[new_location] += 1;
    end

    #7  Resource Growth (R)
    if draw_event >= pr_line[6] && draw_event < pr_line[7]

      #Choose a random resource
      id = rand(collect(1:NR));
      location = Rloc_vec[id];
      nn = ipbc(location,L);
      #What is the LOCAL DENSITY OF RESOURCES?
      nn_id = Rsite_vec[nn];
      LS = sum(nn_id)/length(nn_id); #Local density value

      #if there are any open spots (such that LS > 0)
      if LS > 0
        #Which nn_id == 1? Grab a random one
        open_sites = nn[find(x->x==1,nn_id)];
        new_location = rand(open_sites);

        #place new resource at new location
        push!(Rloc_vec,new_location);
        Rsite_vec[new_location] = 1;

        #Update
        NR = NR + 1;
      end
    end

    #8  Resource consumption (R)
    if draw_event >= pr_line[7] && draw_event < pr_line[8]
      #Draw a random R postion
      id = rand(collect(1:NR));
      location = Rloc_vec[id];

      #Are there ANY F + H in nearest neighbor positions?
      #IF THERE ARE > 0 NF + NH in local neighborhood, then resource is consumed!
      nn = ipbc(location,L);
      nn_idF = Fsite_vec[nn];
      nn_idH = Hsite_vec[nn];
      if sum(nn_idF) + sum(nn_idH) > 0

        #Delete from vector
        # deleteat!(Rind_vec,id);
        deleteat!(Rloc_vec,id);
        Rsite_vec[location] = 0;

        #Update
        NR = NR - 1;
        #Update noresourcesites (add site that resource was eliminated from)
        push!(noresourcesites,location);

      end

    end


    #Recalculate the total # of the foragers and resources
    N = NF + NH + NR;

    #Recalculate the densities of each
    F = NF/S_func;
    H = NH/S_func;
    R = NR/S_func;

    prop = [F,H,R];

    #push!(prop_out,copy(prop));
    prop_out = hcat(prop_out,prop);

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

    if write_out == true
      #Export for Animation
      name = "$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/output/" * "lattice" * @sprintf("%i",tic) * ".csv";
      R_array = Array{Int64}((L+2),(L+2));
      R_array[1:S] = Rsite_vec;
      writedlm(name,R_array,',');
      #R_array_db = DataFrame(R_array);
      #writetable("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/output/name.csv",R_array);
    end

    #ERRORS
    #Break loop if extinction occurs
    if N < 1
      println("Extinction has occured at t=",round(t,2)," and loop ",tic)
      break
    end
    #Break loop if resources go extinct and the other populations run away
    if NR == 0 && (NF + NH) > S
      println("Runaway growth has occured at t=",round(t,2)," and loop ",tic)
      break
    end



  end #end while loop over t

  println("Simulation successful at t= ",round(t,2)," and loop= ",tic)

  # return ind_out,loc_out,time_out,prop_out,N_out
  return time_out,prop_out,N_out
end #end function
