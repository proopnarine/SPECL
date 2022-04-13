function species_level_network(adjacency_matrix,metanetwork,γ,SLN_matrix=nothing,species_statistics=nothing)

    include("./SLN_maker.jl")
    include("./r_no_prey.jl")
    
    P = CSV.read(metanetwork,DataFrame)
    A = readdlm(adjacency_matrix)
    γ = γ

    begin
        # calculate metanetwork diversity
        #no. of guilds
        no_guilds = size(P,1)
        println("\nNo. of guilds, G = ", no_guilds)
        #calculate number of species
        S = sum(P[:,:G])
        no_species = S[1]
        println("No. of species, S = ", no_species)
    end

    begin
        # calculate no. of prey species per guild
        P[:,:no_prey] .= 0.0
        #print(P)
        for i = 1:no_guilds
	    for j = 1:no_guilds
	        if A[i,j] == 1
	            P[i,:no_prey] = P[i,:no_prey] + P[j,:G]
	        end
                if A[j,i] == 1
	            P[i,:no_preds] = P[i,:no_preds] + P[j,:G]
	        end
	    end
        end
        #println("\n", P,"\n")
    end

    # construct guild x species array
    meta_SLN = SLN_maker(A,P,no_guilds,no_species);

    # dataframe for species data
    species = DataFrame(sp_name = Int64[], guild = String[], guild_no = Int64[], g_richness = Int64[], g_no_prey = Int64[], sp_no_prey = Int64[], sp_no_preds = Int64[], sp_long_chain = Int64[], sp_ntp = Float64[])
    
    # push guild data
    begin
        tally1 = [1]
        for i = 1:no_guilds
            guild_richness = P[i,:G]
            for j = 1:guild_richness
                push!(species, [tally1[1], P[i,:Guild], i, P[i,:G], P[i,:no_prey], 0, 0, 0, 0])
                tally1[1] = tally1[1] + 1
            end
        end
    end

    # initial species no. of prey; uses in-degree distribution
    for i = 1:size(species,1)
        species[i,:sp_no_prey] = r_no_prey(species[i,:g_no_prey],γ)
    end
    #print(species)

    # select species-specific prey and generate species A matrix
    sp_A = zeros(Int64,no_species,no_species)
    for i = 1:no_species
        #vector of species prey indices
        guild_prey = Int64[]
        for j = 1:no_species
            if meta_SLN[species[i,:guild_no],j] == 1
                push!(guild_prey,j)
            end
        end
        #randomize prey vector
        shuffle!(guild_prey)
        #pick species prey and predators
        for k = 1:species[i,:sp_no_prey]
            sp_A[i,guild_prey[k]] = 1
        end
    end

    #calculate no. preds, or out-degree
    for i = 1:no_species
        out_degree = 0
        for j = 1:no_species
            if sp_A[j,i] == 1
                out_degree += 1
            end
        end
        species[i,:sp_no_preds] = out_degree
    end

    #initialize pathways matrix
    paths = Array{Int64}(undef,no_species,no_species)
    paths = deepcopy(sp_A)
    #longest possible pathway
    P_max = no_guilds - 1
    #set initial longest path for each species
    for i = 1:no_species
        if species[i,:sp_no_prey] > 0
            species[i,:sp_long_chain] = 1
        end
    end

    #calculate pathways by raising binary adjacency matrix to pathway lengths
    for i = 1:P_max
        A2 = sp_A^i
        for j = 1:no_species
            for k = 1:no_species
                #if path now exists between species
                if paths[j,k]==0 && A2[j,k]>0
                   #update the pathways matrix
                   paths[j,k] = i
                end
            end
            #list as longest chain if one exists
            if sum(paths[j,:]) != 0
                species[j,:sp_long_chain] = maximum(paths[j,:])
            end
        end       
        if sum(A2)==0
            break
        end
    end

    #calculate ntps
    #build vector of primary producers
    prods = Int64[]
    for i = 1:no_species
        if species[i,:sp_no_prey] == 0
            push!(prods,i)
        end
    end
    #calculate path length of prey to producers
    for i = 1:no_species
        #if producer
        if species[i,:sp_no_prey]==0
            species[i,:sp_ntp] = 1
        elseif species[i,:sp_no_prey]>0
            #else if consumer
            #list prey
            its_prey = Int64[]
            path_length = 0
            no_paths = 0
            for j = 1:no_species
                #if species is producer prey of i
                if paths[i,j] == 1 && species[j,:sp_no_prey] == 0
                    no_paths+=1
                end
                #if species is consumer prey of i
                if paths[i,j] == 1 && species[j,:sp_no_prey] > 0
                    #record path lengths to producers
                    for k = 1:no_species
                        if paths[j,k]!=0 && species[k,:sp_no_prey]==0
                            path_length = path_length + paths[j,k]
                            no_paths+=1
                        end
                    end
                end 
            end
            #if herbivore
            if path_length==0
                species[i,:sp_ntp] = 2.0
            elseif path_length > 0
                #if not herbivore
                species[i,:sp_ntp] = 2.0 + (Float64(path_length)/Float64(no_paths))
                #println(species[i,6])
            end
        end
    end
    
    # output
    print(species)
    
    # save output SLN
    if SLN_matrix!=nothing
        writedlm(SLN_matrix,sp_A)
    end

    # save species dataframe
    if species_statistics!=nothing
        CSV.write(species_statistics,species)
    end

    # end function
    return nothing
end


    
