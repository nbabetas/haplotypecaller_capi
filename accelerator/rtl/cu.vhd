library ieee;
  use ieee.std_logic_1164.all;
  use ieee.numeric_std.all;

library work;
  use work.functions.all;
  use work.psl.all;
  use work.cu_package.all;
  use work.dma_package.all;
  use work.pairhmm_package.all;
  use work.pe_package.all;

library xil_defaultlib;

entity cu is
  port (
    i                                       : in  cu_in;
    o                                       : out cu_out
  );
end entity cu;

architecture logic of cu is
  signal q, r                              : cu_int;
  signal re,re2,re3,re4                     : cu_ext;
  signal qs, rs,qs2,rs2,qs3, rs3,qs4,rs4    : cu_sched              := cu_sched_empty;
  signal if_rv                              : unsigned(31 downto 0) := (others => '0');
  signal of_wa                              : unsigned(31 downto 0) := (others => '0');
  signal s1,s2,s3,s4						: std_logic;
  signal delay								: unsigned(3 downto 0) :=(others => '0'); 
  signal temp,temp2							: std_logic_vector(DMA_DATA_WIDTH - 1 downto 0);
  signal read_valid							: std_logic;
begin

---------------------------------------------------------------------------------------------------
--    ____        _       _       _                     _
--   |  _ \      | |     | |     | |                   | |
--   | |_) | __ _| |_ ___| |__   | |     ___   __ _  __| | ___ _ __
--   |  _ < / _` | __/ __| '_ \  | |    / _ \ / _` |/ _` |/ _ \ '__|
--   | |_) | (_| | || (__| | | | | |___| (_) | (_| | (_| |  __/ |
--   |____/ \__,_|\__\___|_| |_| |______\___/ \__,_|\__,_|\___|_|
---------------------------------------------------------------------------------------------------
  loader_comb : process(all)
    variable v                              : cu_int;
    variable t                              : std_logic;
  begin
--------------------------------------------------------------------------------------------------- default assignments 
    v                                       := r;
    v.o.read.valid                          := '0';
    v.o.write.request.valid                 := '0';
    v.o.write.data.valid                    := '0';

    v.read_wren                             := '0';
    v.hapl_wren                             := '0';
    v.p_fifo_en                             := '0';

    v.read_wren_2                           := '0';
    v.hapl_wren_2                           := '0';
    v.p_fifo_en_2                           := '0';

    v.read_wren_3                           := '0';
    v.hapl_wren_3                           := '0';
    v.p_fifo_en_3                           := '0';

    v.read_wren_4                           := '0';
    v.hapl_wren_4                           := '0';
    v.p_fifo_en_4                           := '0';

    v.outfifo_rd                            := '0';
    v.outfifo_rd_2                          := '0';
    v.outfifo_rd_3                          := '0';
    v.outfifo_rd_4                          := '0';
--------------------------------------------------------------------------------------------------- state machine
	 			
				
	
     case r.state is
      -- Idle state, to wait until the WED has been loaded.
      when idle     =>
        -- Reset input FIFO when idle.
        v.fifos_rst                         := '1';

        -- When start signal is received:
        if i.start then
          v.state                           := load_init;
          v.wed                             := i.wed;
          v.resultcachelines                := u(idx(i.wed.batches) * CU_BATCH_RESULT_CACHELINES,32);
          v.o.done                          := '0';
          v.fifos_rst                       := '0';
		  s1							    <= '1';
		  
          -- Get the initial values for the first batch
          read_cacheline(v.o.read, i.wed.source,0);
          v.wed.source                      := i.wed.source + PSL_CACHELINE_SIZE;
        end if;
      
      -- State to register the initial values of the first row of D
      when load_init =>
        if i.read.valid = '1' then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits.x_size                    := u(i.read.data(575 downto 544));
          v.inits.x_padded                  := u(i.read.data(607 downto 576));
          v.inits.y_size                    := u(i.read.data(639 downto 608));
          v.inits.y_padded                  := u(i.read.data(671 downto 640));
          v.inits.x_bppadded                := u(i.read.data(703 downto 672));
          
          v.state                           := request_data;
        end if;
      
      -- Request all data
      when request_data =>
        -- Reset all counters etc...
        v.x_reads                           := (others => '0');
        v.y_reads                           := (others => '0');
        v.p_reads                           := (others => '0');
        v.pair                              := (others => '0');
        v.filled                            := '0';
              
        read_bytes(v.o.read, r.wed.source, v.inits.batch_bytes(DMA_SIZE_WIDTH-1 downto 0),0);
        v.wed.source                        := r.wed.source + v.inits.batch_bytes;
        v.state                           := loadx;
        
        -- Each pair gets a quarter of a cacheline to put the result in
        write_cachelines(v.o.write.request, r.wed.destination, CU_BATCH_RESULT_CACHELINES,0);
        v.wed.destination                 := r.wed.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;

        -- Determine RAM location and flip RAM bit
        v.ram                             := not(r.ram);        
        v.read_addr                       := u(idx(u(v.ram)) * CU_RAM_ADDRS_PER_BATCH,7);
        v.hapl_addr                       := u(idx(u(v.ram)) * CU_RAM_ADDRS_PER_BATCH,7);        

      -- Load the values for the read X
      when loadx =>
        -- Store the bases of the read in the RAM
        if i.read.valid = '1'  then
          v.x_reads                         := r.x_reads + 1;
          v.read_addr                       := r.read_addr + 1;
          v.read_wren                       := '1';
         -- v.read_data                       := i.read.data;
        end if;

        -- If all (padded) bases of all reads are completely loaded,
        -- go to the next state to load the haplotypes
        if v.x_reads = v.inits.x_bppadded / CU_BASES_PER_CACHELINE then
          v.state                           := loady;
        end if;

      -- Load the values for the haplotype Y
      when loady =>
        -- Store the bases of the haplotype in the RAM
        if i.read.valid = '1' then
          v.y_reads                         := r.y_reads + 1;
          v.hapl_addr                       := r.hapl_addr + 1;
          v.hapl_wren                       := '1';
         -- v.hapl_data                       := i.read.data;
        end if;

        -- If all bases of all haplotypes are completely loaded,
        -- go to the next state to load the probabilities
        if v.y_reads = v.inits.y_padded / CU_BASES_PER_CACHELINE then
          v.state                           := streamps;
		   v.wed.batches                     := r.wed.batches - 1;
        end if;

      -- State to stream in the probabilities
      when streamps =>
        -- Keep count of howmany cachelines have come in from the probabilities array.
        v.p_reads                           := r.p_reads + u(i.read.valid);
        -- Enable FIFO read if read is valid:
        v.p_fifo_en                         := i.read.valid;

        -- If we have loaded all the (padded) probabilities of the batch into the FIFO's
        if idx(v.p_reads) = v.inits.x_padded * CU_CACHELINES_PER_PROB then          
         
		  v.state                           := loadnextinit;
          v.p_reads                         := (others => '0');
        end if;        
        -- Copy the initials into the scheduler registers to make room for the new inits.
        v.sched       := v.inits;
        v.sched_array := v.initial_array;

     

      -- Load next batch information
      when loadnextinit => 
        if i.read.valid = '1'  then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits.x_size                    := u(i.read.data(575 downto 544));
          v.inits.x_padded                  := u(i.read.data(607 downto 576));
          v.inits.y_size                    := u(i.read.data(639 downto 608));
          v.inits.y_padded                  := u(i.read.data(671 downto 640));
          v.inits.x_bppadded                := u(i.read.data(703 downto 672));	 
          v.state                           :=launch;
		 
        end if;

		when launch   =>
		v.o.give_2						:='1';
        -- Disable FIFO writes
        v.p_fifo_en      := '0';
        -- If we told the scheduler to start a new batch
        if r.filled = '1' then
          -- And it's not idle anymore, it's busy with the new batch
          if rs.state /= idle then
            -- We can reset the filled bit to 0.
            v.filled                        := '0';
			
            -- If there is still work to do:
            if v.wed.batches /= 0  then 
           	   -- We can start loading a new batch
            	 	 v.state                           := request_data;  
					
			else
             v.state                           := get_wed;
			  			 
			        	 

			 s1	<='0';
            end if;
          end if;
        end if;

		
        -- If the scheduler is idle
        if rs.state = idle and v.state/= done then
          -- We can signal the scheduler to start processing:
          v.filled                          := '1';
        end if;
	 when get_wed =>
			v.wed2						    := i.wed;
			v.batches2						:= i.wed.batches;
			v.state							:= done;	
			
      -- State where we wait for the scheduler to stop
      when done     =>
		v.start_2							  := '1';
        if  rs.state = idle
        and r.resultcachelines = 0
        then			 
         	 write_byte (v.o.write, i.wed_loc1,slv(x"01"),0);
			 v.state                           := done_2;	    
        end if;
	    -- Each pair gets a quarter of a cacheline to put the result in
	  when done_2 =>
		  v.start_2					    := '0';		
		  if i.write.valid='1'  then	   
          v.o.done                          := '1';
		  v.state                           := done_3;	
		  write_cachelines(v.o.write.request, r.wed2.destination, CU_BATCH_RESULT_CACHELINES,0);
          v.wed2.destination                 := r.wed2.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;
		  v.batches2					       := r.batches2 - 1;	 
		end if;	 
	 when done_3 =>
		  if r.batches2 = 0 then		
			v.state 							:=idle;  
		else
			write_cachelines(v.o.write.request, r.wed2.destination, CU_BATCH_RESULT_CACHELINES,0);
            v.wed2.destination                 := r.wed2.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;
			v.batches2					       := r.batches2 - 1;
		end if;
      when others => null;

end case;
---------------------------------------READ SECOND BATCH---------------------------------------------------
 	case r.state_2 is
      when idle     =>
        -- Reset input FIFO when idle.
        v.fifos_rst_2                         := '1';

        -- When start signal is received:
        if r.start_2='1'  then
          v.state_2                           := load_init;          
	      	  
          v.resultcachelines_2              := u(idx(r.wed2.batches) * CU_BATCH_RESULT_CACHELINES,32);
          v.fifos_rst_2                       := '0';
          v.o.done_2                          := '0';
		   s2	<='1';
          -- Get the initial values for the first batch
		   read_cacheline(v.o.read, r.wed2.source,0);
           v.wed2.source                      := r.wed2.source + PSL_CACHELINE_SIZE;
        end if;
      when load_init =>
        if i.read.valid = '1' then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array_2(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits_2.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits_2.x_size                    := u(i.read.data(575 downto 544));
          v.inits_2.x_padded                  := u(i.read.data(607 downto 576));
          v.inits_2.y_size                    := u(i.read.data(639 downto 608));
          v.inits_2.y_padded                  := u(i.read.data(671 downto 640));
          v.inits_2.x_bppadded                := u(i.read.data(703 downto 672));
          
          v.state_2                           :=  request_data;
        end if;
      
      -- Request all data
      when request_data =>
        -- Reset all counters etc...
        v.x_reads_2                           := (others => '0');
        v.y_reads_2                           := (others => '0');
        v.p_reads_2                           := (others => '0');
        v.pair_2                              := (others => '0');
        v.filled_2                            := '0';
              
        read_bytes(v.o.read, r.wed2.source, v.inits_2.batch_bytes(DMA_SIZE_WIDTH-1 downto 0),0);
        v.wed2.source                        := r.wed2.source + v.inits_2.batch_bytes;
         v.state_2                           :=  loadx;
        
     
        
	
        -- Determine RAM location and flip RAM bit
        v.ram_2                             := not(r.ram_2);
        
        v.read_addr_2                       := u(idx(u(v.ram_2)) * CU_RAM_ADDRS_PER_BATCH,7);
        v.hapl_addr_2                       := u(idx(u(v.ram_2)) * CU_RAM_ADDRS_PER_BATCH,7);

	  	
			
		
      -- Load the values for the read X
      when loadx =>
        -- Store the bases of the read in the RAM
        if i.read.valid = '1' then
          v.x_reads_2                         := r.x_reads_2 + 1;
          v.read_addr_2                       := r.read_addr_2 + 1;
          v.read_wren_2                       := '1';
          end if;							  

        -- If all (padded) bases of all reads are completely loaded,
        -- go to the next state to load the haplotypes
        if v.x_reads_2 = v.inits_2.x_bppadded / CU_BASES_PER_CACHELINE then
           v.state_2                           :=  loady;
        end if;

      -- Load the values for the haplotype Y
      when loady =>
        -- Store the bases of the haplotype in the RAM
        if i.read.valid = '1'  then
          v.y_reads_2                         := r.y_reads_2 + 1;
          v.hapl_addr_2                       := r.hapl_addr_2 + 1;
          v.hapl_wren_2                       := '1';
        end if;

        -- If all bases of all haplotypes are completely loaded,
        -- go to the next state to load the probabilities
        if v.y_reads_2 = v.inits_2.y_padded / CU_BASES_PER_CACHELINE then
           v.state_2                           :=  streamps;
		   v.wed2.batches                     := r.wed2.batches - 1;
        end if;

      -- State to stream in the probabilities
      when streamps =>
        -- Keep count of howmany cachelines have come in from the probabilities array.
        v.p_reads_2                           := r.p_reads_2 + u(i.read.valid);

        -- Enable FIFO read if read is valid:
        v.p_fifo_en_2                         := i.read.valid;

        -- If we have loaded all the (padded) probabilities of the batch into the FIFO's
        if idx(v.p_reads_2) = v.inits_2.x_padded * CU_CACHELINES_PER_PROB then
          
		  v.state_2                           := loadnextinit;
          v.p_reads_2                         := (others => '0');         
        end if;
        
        -- Copy the initials into the scheduler registers to make room for the new inits.
        v.sched_2       := v.inits_2;
        v.sched_array_2 := v.initial_array_2;


		  

      -- Load next batch information
      when loadnextinit => 
        if i.read.valid = '1'  then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array_2(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits_2.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits_2.x_size                    := u(i.read.data(575 downto 544));
          v.inits_2.x_padded                  := u(i.read.data(607 downto 576));
          v.inits_2.y_size                    := u(i.read.data(639 downto 608));
          v.inits_2.y_padded                  := u(i.read.data(671 downto 640));
          v.inits_2.x_bppadded                := u(i.read.data(703 downto 672));		
          v.state_2                           :=  launch;
        end if;

		when launch   =>
		 v.o.give_3						:='1';
         v.o.give_2                                           :='0';
        -- Disable FIFO writes
        v.p_fifo_en_2      := '0';		
        -- If we told the scheduler to start a new batch
        if r.filled_2 = '1' then
          -- And it's not idle anymore, it's busy with the new batch
          if rs2.state /= idle then
            -- We can reset the filled bit to 0.
            v.filled_2                        := '0';
            -- If there is still work to do:
            if v.wed2.batches /= 0   then				
              -- We can start loading a new batch
			     	 v.state_2                           :=  request_data; 	
					
								       
            else 				  
             	   v.state_2                           := get_wed;
				 s2	<='0';
            end if;
          end if;
        end if;
		
        -- If the scheduler is idle
        if rs2.state = idle and  v.state_2/= done then
          -- We can signal the scheduler to start processing:
          v.filled_2                          := '1';
        end if;
		when get_wed =>
			v.wed3							:=i.wed;
			v.batches3						:=i.wed.batches;
			v.state_2		  				:=done;
			--v.start_3						  := '1';
      -- State where we wait for the scheduler to stop
      when done     =>
		v.start_3		 := '1';
        if  rs2.state = idle
        and r.resultcachelines_2 = 0
        then
			
		  	write_byte (v.o.write, r.wed.wed2,slv(x"01"),0);	
           v.state_2                           :=  done_2;
		end if;
	  when done_2 =>	
		  v.start_3		 := '0';			
		  if i.write.valid='1'   then			  
          	v.o.done_2                          := '1';
		  	v.state_2                           := done_3;
			write_cachelines(v.o.write.request, r.wed3.destination, CU_BATCH_RESULT_CACHELINES,0);
            v.wed3.destination                 := r.wed3.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;
			v.batches3						   := r.batches3 - 1;
		end if;	 
	 when done_3 =>
		if r.batches3= 0 then
			v.state_2 							:=idle;
		else
			write_cachelines(v.o.write.request, r.wed3.destination, CU_BATCH_RESULT_CACHELINES,0);
            v.wed3.destination                 := r.wed3.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;
			v.batches3						   := r.batches3 - 1;
		 end if;	
		
	when others => null;
end case;
---------------------------------------READ THIRD BATCH---------------------------------------------------
      -- State to register the initial values of the first row of D
 case r.state_3 is
      when idle     =>
        -- Reset input FIFO when idle.
        v.fifos_rst_3                        := '1';

        -- When start signal is received:
        if r.start_3='1'then
          v.state_3                           :=  load_init;          
	     -- v.wed3						    := i.wed3;	  
          v.resultcachelines_3              := u(idx(r.wed3.batches) * CU_BATCH_RESULT_CACHELINES,32);
          v.fifos_rst_3                     := '0';
           v.o.done_3                          := '0';
		   s3	<='1';
          -- Get the initial values for the first batch
			read_cacheline(v.o.read, r.wed3.source,0);
          v.wed3.source                      := r.wed3.source + PSL_CACHELINE_SIZE;	
        end if;
     when load_init =>
        if i.read.valid = '1' then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array_3(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits_3.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits_3.x_size                    := u(i.read.data(575 downto 544));
          v.inits_3.x_padded                  := u(i.read.data(607 downto 576));
          v.inits_3.y_size                    := u(i.read.data(639 downto 608));
          v.inits_3.y_padded                  := u(i.read.data(671 downto 640));
          v.inits_3.x_bppadded                := u(i.read.data(703 downto 672));
          
        v.state_3                           :=  request_data;
        end if;
      
      -- Request all data
      when request_data =>
        -- Reset all counters etc...
        v.x_reads_3                           := (others => '0');
        v.y_reads_3                           := (others => '0');
        v.p_reads_3                           := (others => '0');
        v.pair_3                              := (others => '0');
        v.filled_3                            := '0';
              
        read_bytes(v.o.read, r.wed3.source, v.inits_3.batch_bytes(DMA_SIZE_WIDTH-1 downto 0),0);
        v.wed3.source                        := r.wed3.source + v.inits_3.batch_bytes;
        v.state_3                           :=  loadx;
        
        -- Each pair gets a quarter of a cacheline to put the result in
        
	
        -- Determine RAM location and flip RAM bit
        v.ram_3                             := not(r.ram_3);
        
        v.read_addr_3                      := u(idx(u(v.ram_3)) * CU_RAM_ADDRS_PER_BATCH,7);
        v.hapl_addr_3                       := u(idx(u(v.ram_3)) * CU_RAM_ADDRS_PER_BATCH,7);

        

      -- Load the values for the read X
      when loadx =>
        -- Store the bases of the read in the RAM
        if i.read.valid = '1'  then
          v.x_reads_3                         := r.x_reads_3 + 1;
          v.read_addr_3                       := r.read_addr_3 + 1;
          v.read_wren_3                       := '1';
       --   v.read_data_3                       := i.read.data;
        end if;

        -- If all (padded) bases of all reads are completely loaded,
        -- go to the next state to load the haplotypes
        if v.x_reads_3 = v.inits_3.x_bppadded / CU_BASES_PER_CACHELINE then
          v.state_3                           :=  loady;
        end if;

      -- Load the values for the haplotype Y
      when loady =>
        -- Store the bases of the haplotype in the RAM
        if i.read.valid = '1'  then
          v.y_reads_3                         := r.y_reads_3 + 1;
          v.hapl_addr_3                       := r.hapl_addr_3 + 1;
          v.hapl_wren_3                       := '1';
         -- v.hapl_data_3                       := i.read.data;
        end if;

        -- If all bases of all haplotypes are completely loaded,
        -- go to the next state to load the probabilities
        if v.y_reads_3 = v.inits_3.y_padded / CU_BASES_PER_CACHELINE then
          v.state_3                           :=  streamps;
		  v.wed3.batches                     := r.wed3.batches - 1;
        end if;

      -- State to stream in the probabilities
      when streamps =>
        -- Keep count of howmany cachelines have come in from the probabilities array.
        v.p_reads_3                           := r.p_reads_3 + u(i.read.valid);
        -- Enable FIFO read if read is valid:
        v.p_fifo_en_3                         := i.read.valid;

        -- If we have loaded all the (padded) probabilities of the batch into the FIFO's
        if idx(v.p_reads_3) = v.inits_3.x_padded * CU_CACHELINES_PER_PROB then          
		  v.state_3                           := loadnextinit;
          v.p_reads_3                         := (others => '0');
        end if;
        
        -- Copy the initials into the scheduler registers to make room for the new inits.
        v.sched_3       := v.inits_3;
        v.sched_array_3 := v.initial_array_3;

		  
      -- Load next batch information
      when loadnextinit => 
        if i.read.valid = '1'  then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array_3(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits_3.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits_3.x_size                    := u(i.read.data(575 downto 544));
          v.inits_3.x_padded                  := u(i.read.data(607 downto 576));
          v.inits_3.y_size                    := u(i.read.data(639 downto 608));
          v.inits_3.y_padded                  := u(i.read.data(671 downto 640));
          v.inits_3.x_bppadded                := u(i.read.data(703 downto 672));		
          v.state_3                           :=  launch;
        end if;

		when launch   =>
		v.o.give_4						:='1';
        v.o.give_3                                           :='0';	
        -- Disable FIFO writes
        v.p_fifo_en_3      := '0';		
        -- If we told the scheduler to start a new batch
        if r.filled_3 = '1' then
          -- And it's not idle anymore, it's busy with the new batch
          if rs3.state /= idle then
            -- We can reset the filled bit to 0.
            v.filled_3                        := '0';
            -- If there is still work to do:
            if v.wed3.batches /= 0   then				
              -- We can start loading a new batch
			     	v.state_3                           :=  request_data; 
									       
            else 				  
             	  v.state_3                           :=  get_wed;
				  
				  s3	<='0';
            end if;
          end if;
        end if;

        -- If the scheduler is idle
        if rs3.state = idle and v.state_3 /= done then
          -- We can signal the scheduler to start processing:
          v.filled_3                          := '1';
        end if;
	  when get_wed =>
		
		 v.wed4							:=i.wed;
		  v.batches4						:=i.wed.batches;
		 v.state_3						:= done;
      -- State where we wait for the scheduler to stop
      when done     =>
		v.start_4		 := '1';
        if  rs3.state = idle
        and r.resultcachelines_3 = 0
        then
			
		    write_byte (v.o.write, r.wed2.wed2,slv(x"01"),0);
        	v.state_3                           :=  done_2;
		end if;
	  when done_2 =>
		  v.start_4						  := '0';
		 if i.write.valid='1'  then			  
          	v.o.done_3                          := '1';
			v.state_3                           := done_3;
			write_cachelines(v.o.write.request, r.wed4.destination, CU_BATCH_RESULT_CACHELINES,0 );
            v.wed4.destination                  := r.wed4.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;
			v.batches4						    := r.batches4 - 1;
		 end if;	
	  when done_3 =>
		 if r.batches4 = 0 then
			v.state_3							:=idle;
		 else
			write_cachelines(v.o.write.request, r.wed4.destination, CU_BATCH_RESULT_CACHELINES,0 );
            v.wed4.destination                  := r.wed4.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;
			v.batches4						    := r.batches4 - 1;
		 end if;
	when others => null;
	end case;
---------------------------------------READ FOURTH BATCH---------------------------------------------------
 case r.state_4 is
	  when idle     =>
        -- Reset input FIFO when idle.
        v.fifos_rst_4                       := '1';

        -- When start signal is received:
        if r.start_4='1' then
         v.state_4                           :=  load_init;          
	    --  v.wed4						    := i.wed4;	  
          v.resultcachelines_4              := u(idx(r.wed4.batches) * CU_BATCH_RESULT_CACHELINES,32);
          v.fifos_rst_4                     := '0';
           v.o.done_4                       := '0';
			v.delay		:= (others => '0');
			 s4	<='1';
          -- Get the initial values for the first batch
		  read_cacheline(v.o.read, r.wed4.source,0);
          v.wed4.source                      := i.wed.source + PSL_CACHELINE_SIZE;
        end if;
      when load_init =>
		
        if i.read.valid = '1' then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array_4(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits_4.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits_4.x_size                    := u(i.read.data(575 downto 544));
          v.inits_4.x_padded                  := u(i.read.data(607 downto 576));
          v.inits_4.y_size                    := u(i.read.data(639 downto 608));
          v.inits_4.y_padded                  := u(i.read.data(671 downto 640));
          v.inits_4.x_bppadded                := u(i.read.data(703 downto 672));
          
         v.state_4                           :=request_data;
        end if;
      
      -- Request all data
      when request_data =>
        -- Reset all counters etc...
        v.x_reads_4                           := (others => '0');
        v.y_reads_4                           := (others => '0');
        v.p_reads_4                           := (others => '0');
        v.pair_4                              := (others => '0');
        v.filled_4                            := '0';
              
        read_bytes(v.o.read, r.wed4.source, v.inits_4.batch_bytes(DMA_SIZE_WIDTH-1 downto 0),0);
        v.wed4.source                        := r.wed4.source + v.inits_4.batch_bytes;
        v.state_4                           := loadx;
        
        -- Each pair gets a quarter of a cacheline to put the result in
    --    write_cachelines(v.o.write.request, r.wed4.destination, CU_BATCH_RESULT_CACHELINES,0 );
    --    v.wed4.destination                 := r.wed4.destination + CU_BATCH_RESULT_CACHELINES * PSL_CACHELINE_SIZE;
	
        -- Determine RAM location and flip RAM bit
        v.ram_4                             := not(r.ram_4);
        
        v.read_addr_4                       := u(idx(u(v.ram_4)) * CU_RAM_ADDRS_PER_BATCH,7);
        v.hapl_addr_4                       := u(idx(u(v.ram_4)) * CU_RAM_ADDRS_PER_BATCH,7);

        

      -- Load the values for the read X
      when loadx =>
		
        -- Store the bases of the read in the RAM
        if i.read.valid = '1'  then
          v.x_reads_4                         := r.x_reads_4 + 1;
          v.read_addr_4                       := r.read_addr_4 + 1;
          v.read_wren_4                       := '1';
          --v.read_data_4                       := i.read.data;
        end if;

        -- If all (padded) bases of all reads are completely loaded,
        -- go to the next state to load the haplotypes
        if v.x_reads_4 = v.inits_4.x_bppadded / CU_BASES_PER_CACHELINE then
          v.state_4                           := loady;
        end if;

      -- Load the values for the haplotype Y
      when loady =>
        -- Store the bases of the haplotype in the RAM
        if i.read.valid = '1'   then
          v.y_reads_4                         := r.y_reads_4 + 1;
          v.hapl_addr_4                       := r.hapl_addr_4 + 1;
          v.hapl_wren_4                       := '1';
      --    v.hapl_data_4                       := i.read.data;
        end if;

        -- If all bases of all haplotypes are completely loaded,
        -- go to the next state to load the probabilities
        if v.y_reads_4 = v.inits_4.y_padded / CU_BASES_PER_CACHELINE then
         v.state_4                           :=streamps;
		 v.wed4.batches                      := r.wed4.batches - 1;
        end if;

      -- State to stream in the probabilities
      when streamps =>
        -- Keep count of howmany cachelines have come in from the probabilities array.
        v.p_reads_4                           := r.p_reads_4 + u(i.read.valid);

        -- Enable FIFO read if read is valid:
        v.p_fifo_en_4                         := i.read.valid;

        -- If we have loaded all the (padded) probabilities of the batch into the FIFO's
        if idx(v.p_reads_4) = v.inits_4.x_padded * CU_CACHELINES_PER_PROB then
          
		  v.state_4                           := loadnextinit;
          v.p_reads_4                         := (others => '0');
        end if;
        
        -- Copy the initials into the scheduler registers to make room for the new inits.
        v.sched_4       := v.inits_4;
        v.sched_array_4 := v.initial_array_4;
		  
   
      -- Load next batch information
      when loadnextinit=> 
        if i.read.valid = '1'  then
          for K in 0 to PE_DEPTH - 1 loop
            v.initial_array_4(K)              := i.read.data((K+1) * PE_DW -1 downto K * PE_DW);
          end loop;
          
          v.inits_4.batch_bytes               := u(i.read.data(543 downto 512));
          v.inits_4.x_size                    := u(i.read.data(575 downto 544));
          v.inits_4.x_padded                  := u(i.read.data(607 downto 576));
          v.inits_4.y_size                    := u(i.read.data(639 downto 608));
          v.inits_4.y_padded                  := u(i.read.data(671 downto 640));
          v.inits_4.x_bppadded                := u(i.read.data(703 downto 672));
          v.state_4                           := launch;
		  v.start							  := '1';
        end if;

 -- A new batch is ready to be started
		when launch   =>
        -- Disable FIFO writes
        v.p_fifo_en_4      := '0'; 
        -- If we told the scheduler to start a new batch
        if r.filled_4 = '1' then
          -- And it's not idle anymore, it's busy with the new batch
          if rs4.state /= idle then
            -- We can reset the filled bit to 0.
            v.filled_4                        := '0';
            -- If there is still work to do:
            if v.wed4.batches /= 0
			then
              -- We can start loading a new batch
					v.start_4 						:= '0';		
              		v.state_4                           :=request_data;		
            else 				
             	v.state_4                           := done;	
					 s4	<='0';		
            end if;
          end if;
        end if;

        -- If the scheduler is idle
        if rs4.state = idle and v.state_4/= done then
          -- We can signal the scheduler to start processing:
          v.filled_4                          := '1';
        end if;

      -- State where we wait for the scheduler to stop
      when done     =>
        if  rs4.state = idle
        and r.resultcachelines_4 = 0 
        then	
			write_byte (v.o.write,  r.wed3.wed2,slv(x"01"),0);
			v.state_4                          := done_2;
        end if;
	  when done_2 	=>
		if i.write.valid='1'   then
		v.state_4                           := idle;
        v.o.done_4                          := '1';
		end if;		
      when others => null;

    end case;

--------------------------------------------------------------------------------------------------- output fifo reads

	if      (re.outfifo.c.valid = '1')  -- if there is a valid output
        and (i.write.full(0) = '0')    -- and the write buffer is not full
    then
      --v.outfifo_rd                          := '1';
      write_data                            (v.o.write.data, re.outdata,0);
	end if;	
	if      (re2.outfifo.c.valid = '1')   -- if there is a valid output
        and (i.write.full(0) = '0') and (r.o.done ='1')  -- and the write buffer is not full
    then
     -- v.outfifo_rd_2                          := '1';
      write_data                            (v.o.write.data, re2.outdata,0);
	end if;
	if      (re3.outfifo.c.valid = '1') -- if there is a valid output
        and (i.write.full(0) = '0') and (r.o.done_2 ='1')     -- and the write buffer is not full
    then
    --  v.outfifo_rd_3                          := '1';
      write_data                            (v.o.write.data, re3.outdata,0);
	end if;
	if      (re4.outfifo.c.valid = '1') -- if there is a valid output
        and (i.write.full(0) = '0') and (r.o.done_3 ='1')   -- and the write buffer is not full
    then
     -- v.outfifo_rd_4                          := '1';
      write_data                            (v.o.write.data, re4.outdata,0);
    end if;
    
	
	
	

    -- Keep track of howmany cachelines have been written back.
    if i.write.valid = '1' and r.state /= idle  then
      v.resultcachelines                    := r.resultcachelines - 1;
    
    end if;

    if i.write.valid = '1' and r.state_2 /= idle and r.state =idle then
      v.resultcachelines_2                    := r.resultcachelines_2 - 1;
       
    end if;

	if i.write.valid = '1' and r.state_3 /= idle and r.state_2 = idle then
      v.resultcachelines_3                    := r.resultcachelines_3 - 1;
    
    end if;

    if i.write.valid = '1' and r.state_4 /= idle and r.state_3 = idle then
      v.resultcachelines_4                    := r.resultcachelines_4 - 1;
       
    end if;
--------------------------------------------------------------------------------------------------- status register updates
    if r.o.mmio_regs.fifo = X"0000000000000000" then
      v.o.mmio_regs.fifo(0)                 := re2.outfifo.c.underflow;
      v.o.mmio_regs.fifo(1)                 := re2.outfifo.c.overflow;
    end if;

--------------------------------------------------------------------------------------------------- signals with latency
    -- These signals are either to match latencies of other components or  to cut critical paths up
    v.read_addr1                               := r.read_addr;
    v.hapl_addr1                               := r.hapl_addr;

    v.read_addr1_2                             := r.read_addr_2;
    v.hapl_addr1_2                             := r.hapl_addr_2;

    v.read_addr1_3                             := r.read_addr_3;
    v.hapl_addr1_3                             := r.hapl_addr_3;

    v.read_addr1_4                             := r.read_addr_4;
    v.hapl_addr1_4                             := r.hapl_addr_4;

    if s1 ='1' then 
       v.p_fifodata(0)                         := i.read.data;
       
    end if;
    v.p_fifodata(1)                         := r.p_fifodata(0);


    if s2 ='1' then
       v.p_fifodata_2(0)                       := i.read.data;
      
    end if;
     v.p_fifodata_2(1)                       := r.p_fifodata_2(0);

    if s3='1' then 
       v.p_fifodata_3(0)                         := i.read.data;
       
    end if;
    v.p_fifodata_3(1)                         := r.p_fifodata_3(0);


    if s4='1' then
       v.p_fifodata_4(0)                       := i.read.data;
      
    end if;
     v.p_fifodata_4(1)                       := r.p_fifodata_4(0);

    v.p_fifo_en1                               := r.p_fifo_en;
    v.p_fifo_en1_2                             := r.p_fifo_en_2;
	v.p_fifo_en1_3                             := r.p_fifo_en_3;
    v.p_fifo_en1_4                             := r.p_fifo_en_4;


--------------------------------------------------------------------------------------------------- outputs
    -- drive input registers
    q                       <= v;
    -- outputs
    o                       <= r.o;
  end process;

--------------------------------------------------------------------------------------------------- registers
  loader_reg : process(i.cr)
  begin
    if rising_edge(i.cr.clk) then
      if i.cr.rst then
        cu_reset(r);
      else
        	r                  		<= q;			
			
      end if;
    end if;
  end process;

loader_reg_2 : process(i.cr)
  begin     	
    if rising_edge(i.cr.clk) then
      temp                                    <=i.read.data;
	--  read_valid							  <=i.read.valid;
    end if;
  end process;

--------------------------------------------------------------------------------------------------- Clock generator

  -- In case the kernel has to run slower due to timing constraints not being met, use this to lower the clock frequency
  kernel_clock_gen : entity xil_defaultlib.psl_to_kernel port map (
    clk_psl            => i.cr.clk,
    clk_kernel         => re.clk_kernel
  ); -- keep one clock for all kernel in order to meet all kerners timing constraints

  -- Use this to keep everything in the same clock domain:
  --re.clk_kernel <= i.cr.clk;
	
---------------------------------------------------------------------------------------------------
--    ____                        _____            __  __
--   |  _ \                      |  __ \     /\   |  \/  |
--   | |_) | __ _ ___  ___  ___  | |__) |   /  \  | \  / |___
--   |  _ < / _` / __|/ _ \/ __| |  _  /   / /\ \ | |\/| / __|
--   | |_) | (_| \__ \  __/\__ \ | | \ \  / ____ \| |  | \__ \
--   |____/ \__,_|___/\___||___/ |_|  \_\/_/    \_\_|  |_|___/
---------------------------------------------------------------------------------------------------
  
  -- Connect clocks
  re.haplram.clka <= i.cr.clk;
  re.readram.clka <= i.cr.clk;
  re.haplram.clkb <= re.clk_kernel;
  re.readram.clkb <= re.clk_kernel;

  -- Connect clocks 2
  re2.haplram.clka <= i.cr.clk;
  re2.readram.clka <= i.cr.clk;
  re2.haplram.clkb <= re.clk_kernel;
  re2.readram.clkb <= re.clk_kernel;
  
  -- Connect clocks 3
  re3.haplram.clka <= i.cr.clk;
  re3.readram.clka <= i.cr.clk;
  re3.haplram.clkb <= re.clk_kernel;
  re3.readram.clkb <= re.clk_kernel;

  -- Connect clocks 4
  re4.haplram.clka <= i.cr.clk;
  re4.readram.clka <= i.cr.clk;
  re4.haplram.clkb <= re.clk_kernel;
  re4.readram.clkb <= re.clk_kernel;


------------AYTO KANE META _____
  --temp2 <= temp;
-------------------------------
  -- Convert data to basepair type
  convert_hapl: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re.haplram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;

  convert_read: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re.readram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;
  
  -- Conver data to basepair type 2
  convert_hapl_2: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re2.haplram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;

  convert_read_2: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re2.readram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;

  -- Convert data to basepair type 3
  convert_hapl_3: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re3.haplram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;

  convert_read_3: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re3.readram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;
  
  -- Conver data to basepair type 4
  convert_hapl_4: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re4.haplram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;

  convert_read_4: for K in 0 to PAIRHMM_BASEPAIRS_PER_CACHELINE-1 generate
    re4.readram.dina((K+1) * BP_SIZE - 1 downto K * BP_SIZE) <= slv8bpslv3(temp((K + 1) * 8 - 1 downto (K * 8)));
  end generate;

  -- Connect write addresses
  re.haplram.addra  <= r.hapl_addr1;
  re.readram.addra  <= r.read_addr1;

  -- Connect write addresses 2
  re2.haplram.addra  <= r.hapl_addr1_2;
  re2.readram.addra  <= r.read_addr1_2;

  -- Connect write addresses 3
  re3.haplram.addra  <= r.hapl_addr1_3;
  re3.readram.addra  <= r.read_addr1_3;

  -- Connect write addresses 4
  re4.haplram.addra  <= r.hapl_addr1_4;
  re4.readram.addra  <= r.read_addr1_4;

  -- Connect write enables
  re.haplram.wea    <= r.hapl_wren;
  re.readram.wea    <= r.read_wren;

  -- Connect write enables 2
  re2.haplram.wea    <= r.hapl_wren_2;
  re2.readram.wea    <= r.read_wren_2;

  -- Connect write enables 3
  re3.haplram.wea    <= r.hapl_wren_3;
  re3.readram.wea    <= r.read_wren_3;

  -- Connect write enables 4
  re4.haplram.wea    <= r.hapl_wren_4;
  re4.readram.wea    <= r.read_wren_4;

  -- Connect addresses to read from haplotype
  re.haplram.addrb  <= u(idx(rs.element) + idx(rs.supercolumn) * PAIRHMM_NUM_PES + PAIRHMM_MAX_SIZE * idx(u(rs.ram)), 10);

  -- Connect addresses to read from haplotype 2
  re2.haplram.addrb  <= u(idx(rs2.element) + idx(rs2.supercolumn) * PAIRHMM_NUM_PES + PAIRHMM_MAX_SIZE * idx(u(rs2.ram)), 10);

 -- Connect addresses to read from haplotype 3
  re3.haplram.addrb  <= u(idx(rs3.element) + idx(rs3.supercolumn) * PAIRHMM_NUM_PES + PAIRHMM_MAX_SIZE * idx(u(rs3.ram)), 10);

  -- Connect addresses to read from haplotype 4
  re4.haplram.addrb  <= u(idx(rs4.element) + idx(rs4.supercolumn) * PAIRHMM_NUM_PES + PAIRHMM_MAX_SIZE * idx(u(rs4.ram)), 10);


  -- Connect addresses to read from read
  re.readram.addrb  <= u(idx(rs.basepair) + PAIRHMM_MAX_SIZE * idx(u(rs.ram)), 10);

  -- Connect addresses to read from read 2
  re2.readram.addrb  <= u(idx(rs2.basepair) + PAIRHMM_MAX_SIZE * idx(u(rs2.ram)), 10);

  -- Connect addresses to read from read 3
  re3.readram.addrb  <= u(idx(rs3.basepair) + PAIRHMM_MAX_SIZE * idx(u(rs3.ram)), 10);

  -- Connect addresses to read from read 4
  re4.readram.addrb  <= u(idx(rs4.basepair) + PAIRHMM_MAX_SIZE * idx(u(rs4.ram)), 10);


  -- Haplotype RAM
  hapl_ram   : asym_ram port map (
    clka    =>     re.haplram.clka,
    wea     =>     re.haplram.wea,
    addra   =>     re.haplram.addra,
    dina    =>     re.haplram.dina,
    clkb    =>     re.haplram.clkb,
    addrb   =>     re.haplram.addrb,
    doutb   =>     re.haplram.doutb
  );

    -- Haplotype RAM 2
  hapl_ram_2   : asym_ram port map (
    clka    =>     re2.haplram.clka,
    wea     =>     re2.haplram.wea,
    addra   =>     re2.haplram.addra,
    dina    =>     re2.haplram.dina,
    clkb    =>     re2.haplram.clkb,
    addrb   =>     re2.haplram.addrb,
    doutb   =>     re2.haplram.doutb
  );

  -- Haplotype RAM 3
  hapl_ram_3   : asym_ram port map (
    clka    =>     re3.haplram.clka,
    wea     =>     re3.haplram.wea,
    addra   =>     re3.haplram.addra,
    dina    =>     re3.haplram.dina,
    clkb    =>     re3.haplram.clkb,
    addrb   =>     re3.haplram.addrb,
    doutb   =>     re3.haplram.doutb
  );

    -- Haplotype RAM 4
  hapl_ram_4   : asym_ram port map (
    clka    =>     re4.haplram.clka,
    wea     =>     re4.haplram.wea,
    addra   =>     re4.haplram.addra,
    dina    =>     re4.haplram.dina,
    clkb    =>     re4.haplram.clkb,
    addrb   =>     re4.haplram.addrb,
    doutb   =>     re4.haplram.doutb
  );
  -- Read RAM
  read_ram  : asym_ram port map (
    clka    =>     re.readram.clka,
    wea     =>     re.readram.wea,
    addra   =>     re.readram.addra,
    dina    =>     re.readram.dina,
    clkb    =>     re.readram.clkb,
    addrb   =>     re.readram.addrb,
    doutb   =>     re.readram.doutb
  );

 -- Read RAM 2
  read_ram_2  : asym_ram port map (
    clka    =>     re2.readram.clka,
    wea     =>     re2.readram.wea,
    addra   =>     re2.readram.addra,
    dina    =>     re2.readram.dina,
    clkb    =>     re2.readram.clkb,
    addrb   =>     re2.readram.addrb,
    doutb   =>     re2.readram.doutb
  );

  -- Read RAM 3
  read_ram_3  : asym_ram port map (
    clka    =>     re3.readram.clka,
    wea     =>     re3.readram.wea,
    addra   =>     re3.readram.addra,
    dina    =>     re3.readram.dina,
    clkb    =>     re3.readram.clkb,
    addrb   =>     re3.readram.addrb,
    doutb   =>     re3.readram.doutb
  );

 -- Read RAM 4
  read_ram_4  : asym_ram port map (
    clka    =>     re4.readram.clka,
    wea     =>     re4.readram.wea,
    addra   =>     re4.readram.addra,
    dina    =>     re4.readram.dina,
    clkb    =>     re4.readram.clkb,
    addrb   =>     re4.readram.addrb,
    doutb   =>     re4.readram.doutb
  );

---------------------------------------------------------------------------------------------------
--    ______ _____ ______ ____
--   |  ____|_   _|  ____/ __ \
--   | |__    | | | |__ | |  | |___
--   |  __|   | | |  __|| |  | / __|
--   | |     _| |_| |   | |__| \__ \
--   |_|    |_____|_|    \____/|___/
----------------------------------------------------------------------------------------------  -----

--------------------------------------------------------------------------------------------------- Input FIFO

    -- Connect reset from the pair loader
    re.pfifo.c.rst     <= r.fifos_rst;
    re2.pfifo.c.rst     <= r.fifos_rst_2;
	re3.pfifo.c.rst     <= r.fifos_rst_3;
    re4.pfifo.c.rst     <= r.fifos_rst_4;

    -- Connect read enables from the scheduler
    re.pfifo.c.rd_en   <= rs.fifo_rd_en;
    re2.pfifo.c.rd_en  <= rs2.fifo_rd_en;
    re3.pfifo.c.rd_en  <= rs3.fifo_rd_en;
    re4.pfifo.c.rd_en  <= rs4.fifo_rd_en;	 
    -- Connect FIFO inputs:
    re.pfifo.c.wr_en    <= r.p_fifo_en1;
    re2.pfifo.c.wr_en   <= r.p_fifo_en1_2;
	re3.pfifo.c.wr_en   <= r.p_fifo_en1_3;
    re4.pfifo.c.wr_en   <= r.p_fifo_en1_4;

    -- Swap the first floats to go into the fifo first
    floatswap: for K in 0 to DMA_DATA_WIDTH / PE_DW - 1 generate
      re.pfifo.din((K+1)*PE_DW-1 downto K*PE_DW) <=  r.p_fifodata(1)(
                                                          DMA_DATA_WIDTH -  K    * PE_DW - 1
                                                          downto
                                                          DMA_DATA_WIDTH - (K+1) * PE_DW
                                                        );

    end generate;

    floatswap_2: for K in 0 to DMA_DATA_WIDTH / PE_DW - 1 generate
    re2.pfifo.din((K+1)*PE_DW-1 downto K*PE_DW) <=  r.p_fifodata_2(1)(
                                                          DMA_DATA_WIDTH -  K    * PE_DW - 1
                                                          downto
                                                          DMA_DATA_WIDTH - (K+1) * PE_DW
                                                        );

    end generate;

    floatswap_3: for K in 0 to DMA_DATA_WIDTH / PE_DW - 1 generate
      re3.pfifo.din((K+1)*PE_DW-1 downto K*PE_DW) <=  r.p_fifodata_3(1)(
                                                          DMA_DATA_WIDTH -  K    * PE_DW - 1
                                                          downto
                                                          DMA_DATA_WIDTH - (K+1) * PE_DW
                                                        );

    end generate;

    -- Swap the first floats to go into the fifo first
    floatswap_4: for K in 0 to DMA_DATA_WIDTH / PE_DW - 1 generate
    re4.pfifo.din((K+1)*PE_DW-1 downto K*PE_DW) <=  r.p_fifodata_4(1)(
                                                          DMA_DATA_WIDTH -  K    * PE_DW - 1
                                                          downto
                                                          DMA_DATA_WIDTH - (K+1) * PE_DW
                                                        );

    end generate;

    -- Connect PAIRHMM inputs 1
    re.pairhmm_in.en                  <= '1';
    re.pairhmm_in.valid               <= rs.valid;
    re.pairhmm_in.cell                <= rs.cell;

    -- Connect PAIRHMM inputs 2
    re2.pairhmm_in.en                  <= '1';
    re2.pairhmm_in.valid               <= rs2.valid;
    re2.pairhmm_in.cell                <= rs2.cell;

    -- Connect PAIRHMM inputs 3
    re3.pairhmm_in.en                  <= '1';
    re3.pairhmm_in.valid               <= rs3.valid;
    re3.pairhmm_in.cell                <= rs3.cell;

    -- Connect PAIRHMM inputs 4
    re4.pairhmm_in.en                  <= '1';
    re4.pairhmm_in.valid               <= rs4.valid;
    re4.pairhmm_in.cell                <= rs4.cell;

    -- Set top left input to 1.0 when this is the first cycle of this pair.

    -- Initial input for first PE
    re.pairhmm_in.mids.dtl             <= rs.initial_array(idx(rs.schedule));
    re2.pairhmm_in.mids.dtl            <= rs2.initial_array(idx(rs2.schedule));
    re3.pairhmm_in.mids.dtl            <= rs3.initial_array(idx(rs3.schedule));
    re4.pairhmm_in.mids.dtl            <= rs4.initial_array(idx(rs4.schedule));
    
    -- Select initial value to travel with systolic array
    re.pairhmm_in.initial             <= rs.initial_array(idx(rs.schedule));

    re.pairhmm_in.mids.itl            <= (others => '0');
    re.pairhmm_in.mids.mtl            <= (others => '0');
    re.pairhmm_in.mids.ml             <= (others => '0');
    re.pairhmm_in.mids.il             <= (others => '0');
    re.pairhmm_in.mids.dl             <= (others => '0');
    re.pairhmm_in.mids.mt             <= (others => '0');
    re.pairhmm_in.mids.it             <= (others => '0');
    re.pairhmm_in.mids.dt             <= (others => '0');

    re.pairhmm_in.emis.distm_simi     <= re.pfifo.dout( 31 downto   0);
    re.pairhmm_in.emis.distm_diff     <= re.pfifo.dout( 63 downto  32);
    re.pairhmm_in.tmis.alpha          <= re.pfifo.dout( 95 downto  64);
    re.pairhmm_in.tmis.beta           <= re.pfifo.dout(127 downto  96);
    re.pairhmm_in.tmis.delta          <= re.pfifo.dout(159 downto 128);
    re.pairhmm_in.tmis.epsilon        <= re.pfifo.dout(191 downto 160);
    re.pairhmm_in.tmis.zeta           <= re.pfifo.dout(223 downto 192);
    re.pairhmm_in.tmis.eta            <= re.pfifo.dout(255 downto 224);

---------------------------------------------------------------------------------------------------2
    re2.pairhmm_in.initial             <= rs2.initial_array(idx(rs2.schedule));

    re2.pairhmm_in.mids.itl            <= (others => '0');
    re2.pairhmm_in.mids.mtl            <= (others => '0');
    re2.pairhmm_in.mids.ml             <= (others => '0');
    re2.pairhmm_in.mids.il             <= (others => '0');
    re2.pairhmm_in.mids.dl             <= (others => '0');
    re2.pairhmm_in.mids.mt             <= (others => '0');
    re2.pairhmm_in.mids.it             <= (others => '0');
    re2.pairhmm_in.mids.dt             <= (others => '0');

    re2.pairhmm_in.emis.distm_simi     <= re2.pfifo.dout( 31 downto   0);
    re2.pairhmm_in.emis.distm_diff     <= re2.pfifo.dout( 63 downto  32);
    re2.pairhmm_in.tmis.alpha          <= re2.pfifo.dout( 95 downto  64);
    re2.pairhmm_in.tmis.beta           <= re2.pfifo.dout(127 downto  96);
    re2.pairhmm_in.tmis.delta          <= re2.pfifo.dout(159 downto 128);
    re2.pairhmm_in.tmis.epsilon        <= re2.pfifo.dout(191 downto 160);
    re2.pairhmm_in.tmis.zeta           <= re2.pfifo.dout(223 downto 192);
    re2.pairhmm_in.tmis.eta            <= re2.pfifo.dout(255 downto 224);
---------------------------------------------------------------------------------------------------3
    re3.pairhmm_in.initial             <= rs3.initial_array(idx(rs3.schedule));

    re3.pairhmm_in.mids.itl            <= (others => '0');
    re3.pairhmm_in.mids.mtl            <= (others => '0');
    re3.pairhmm_in.mids.ml             <= (others => '0');
    re3.pairhmm_in.mids.il             <= (others => '0');
    re3.pairhmm_in.mids.dl             <= (others => '0');
    re3.pairhmm_in.mids.mt             <= (others => '0');
    re3.pairhmm_in.mids.it             <= (others => '0');
    re3.pairhmm_in.mids.dt             <= (others => '0');

    re3.pairhmm_in.emis.distm_simi     <= re3.pfifo.dout( 31 downto   0);
    re3.pairhmm_in.emis.distm_diff     <= re3.pfifo.dout( 63 downto  32);
    re3.pairhmm_in.tmis.alpha          <= re3.pfifo.dout( 95 downto  64);
    re3.pairhmm_in.tmis.beta           <= re3.pfifo.dout(127 downto  96);
    re3.pairhmm_in.tmis.delta          <= re3.pfifo.dout(159 downto 128);
    re3.pairhmm_in.tmis.epsilon        <= re3.pfifo.dout(191 downto 160);
    re3.pairhmm_in.tmis.zeta           <= re3.pfifo.dout(223 downto 192);
    re3.pairhmm_in.tmis.eta            <= re3.pfifo.dout(255 downto 224);
---------------------------------------------------------------------------------------------------4
    re4.pairhmm_in.initial             <= rs4.initial_array(idx(rs4.schedule));

    re4.pairhmm_in.mids.itl            <= (others => '0');
    re4.pairhmm_in.mids.mtl            <= (others => '0');
    re4.pairhmm_in.mids.ml             <= (others => '0');
    re4.pairhmm_in.mids.il             <= (others => '0');
    re4.pairhmm_in.mids.dl             <= (others => '0');
    re4.pairhmm_in.mids.mt             <= (others => '0');
    re4.pairhmm_in.mids.it             <= (others => '0');
    re4.pairhmm_in.mids.dt             <= (others => '0');

    re4.pairhmm_in.emis.distm_simi     <= re4.pfifo.dout( 31 downto   0);
    re4.pairhmm_in.emis.distm_diff     <= re4.pfifo.dout( 63 downto  32);
    re4.pairhmm_in.tmis.alpha          <= re4.pfifo.dout( 95 downto  64);
    re4.pairhmm_in.tmis.beta           <= re4.pfifo.dout(127 downto  96);
    re4.pairhmm_in.tmis.delta          <= re4.pfifo.dout(159 downto 128);
    re4.pairhmm_in.tmis.epsilon        <= re4.pfifo.dout(191 downto 160);
    re4.pairhmm_in.tmis.zeta           <= re4.pfifo.dout(223 downto 192);
    re4.pairhmm_in.tmis.eta            <= re4.pfifo.dout(255 downto 224);

  -- First Word Fall Through FIFO's which can potentially cross clockdomains to a lower clock:
    p_fifo : entity xil_defaultlib.probabilities_fifo port map (
      wr_clk    => i.cr.clk,
      rst       => re.pfifo.c.rst,
      rd_clk    => re.clk_kernel,
      din       => re.pfifo.din,
      wr_en     => re.pfifo.c.wr_en,
      wr_ack    => re.pfifo.c.wr_ack,
      rd_en     => re.pfifo.c.rd_en,
      valid     => re.pfifo.c.valid,
      dout      => re.pfifo.dout,
      full      => re.pfifo.c.full,
      empty     => re.pfifo.c.empty,
      overflow  => re.pfifo.c.overflow,
      underflow => re.pfifo.c.underflow
    );

 p_fifo_2 : entity xil_defaultlib.probabilities_fifo port map (
      wr_clk    => i.cr.clk,
      rst       => re2.pfifo.c.rst,
      rd_clk    => re.clk_kernel,
      din       => re2.pfifo.din,
      wr_en     => re2.pfifo.c.wr_en,
      wr_ack    => re2.pfifo.c.wr_ack,
      rd_en     => re2.pfifo.c.rd_en,
      valid     => re2.pfifo.c.valid,
      dout      => re2.pfifo.dout,
      full      => re2.pfifo.c.full,
      empty     => re2.pfifo.c.empty,
      overflow  => re2.pfifo.c.overflow,
      underflow => re2.pfifo.c.underflow
    );

  -- First Word Fall Through FIFO's which can potentially cross clockdomains to a lower clock:
    p_fifo_3 : entity xil_defaultlib.probabilities_fifo port map (
      wr_clk    => i.cr.clk,
      rst       => re3.pfifo.c.rst,
      rd_clk    => re.clk_kernel,
      din       => re3.pfifo.din,
      wr_en     => re3.pfifo.c.wr_en,
      wr_ack    => re3.pfifo.c.wr_ack,
      rd_en     => re3.pfifo.c.rd_en,
      valid     => re3.pfifo.c.valid,
      dout      => re3.pfifo.dout,
      full      => re3.pfifo.c.full,
      empty     => re3.pfifo.c.empty,
      overflow  => re3.pfifo.c.overflow,
      underflow => re3.pfifo.c.underflow
    );

 p_fifo_4 : entity xil_defaultlib.probabilities_fifo port map (
      wr_clk    => i.cr.clk,
      rst       => re4.pfifo.c.rst,
      rd_clk    => re.clk_kernel,
      din       => re4.pfifo.din,
      wr_en     => re4.pfifo.c.wr_en,
      wr_ack    => re4.pfifo.c.wr_ack,
      rd_en     => re4.pfifo.c.rd_en,
      valid     => re4.pfifo.c.valid,
      dout      => re4.pfifo.dout,
      full      => re4.pfifo.c.full,
      empty     => re4.pfifo.c.empty,
      overflow  => re4.pfifo.c.overflow,
      underflow => re4.pfifo.c.underflow
    );
--------------------------------------------------------------------------------------------------- Output FIFO

  -- Output data that is written back to the memory, goes into the fifo first
  re.outfifo.din  <= ( 31 downto   0 => re.pairhmm.o.score,
                       63 downto  32 => X"00000000", --re.pairhmm.o.score,
                       95 downto  64 => X"00000000", --re.pairhmm.o.score,
                      127 downto  96 => X"00000000"); --re.pairhmm.o.score);

-- Output data that is written back to the memory, goes into the fifo first  SECOND
  re2.outfifo.din  <= ( 31 downto   0 => re2.pairhmm.o.score,
                       63 downto  32 => X"00000000", --re.pairhmm.o.score,
                       95 downto  64 => X"00000000", --re.pairhmm.o.score,
                      127 downto  96 => X"00000000"); --re.pairhmm.o.score);
  -- Output data that is written back to the memory, goes into the fifo first
  re3.outfifo.din  <= ( 31 downto   0 => re3.pairhmm.o.score,
                       63 downto  32 => X"00000000", --re.pairhmm.o.score,
                       95 downto  64 => X"00000000", --re.pairhmm.o.score,
                      127 downto  96 => X"00000000"); --re.pairhmm.o.score);

-- Output data that is written back to the memory, goes into the fifo first  SECOND
  re4.outfifo.din  <= ( 31 downto   0 => re4.pairhmm.o.score,
                       63 downto  32 => X"00000000", --re.pairhmm.o.score,
                       95 downto  64 => X"00000000", --re.pairhmm.o.score,
                      127 downto  96 => X"00000000"); --re.pairhmm.o.score);

  re.outfifo.c.wr_en <= re.pairhmm.o.score_valid;--'1' when (re.pairhmm.o.last.valid='1') and (re.pairhmm.o.last.cell = PE_LAST) else '0';
  re2.outfifo.c.wr_en <= re2.pairhmm.o.score_valid;--'1' when (re.pairhmm.o.last.valid='1') and (re.pairhmm.o.last.cell = PE_LAST) else '0';
  re3.outfifo.c.wr_en <= re3.pairhmm.o.score_valid;--'1' when (re.pairhmm.o.last.valid='1') and (re.pairhmm.o.last.cell = PE_LAST) else '0';
  re4.outfifo.c.wr_en <= re4.pairhmm.o.score_valid;--'1' when (re.pairhmm.o.last.valid='1') and (re.pairhmm.o.last.cell = PE_LAST) else '0';
  -- Is this ugly? Maybe...
  --re.outfifo.c.rd_en <= r.outfifo_rd;
  --re2.outfifo.c.rd_en <= r.outfifo_rd_2; 
  --re3.outfifo.c.rd_en <= r.outfifo_rd_3;
  --re4.outfifo.c.rd_en <= r.outfifo_rd_4; 

  -- This fifo can potentially cross clock domains to the streaming framework.
  outfifo: entity xil_defaultlib.kernel_to_streaming_fifo
    port map (
      wr_clk    => re.clk_kernel,
      rd_clk    => i.cr.clk,
      din       => re.outfifo.din,
      wr_en     => re.outfifo.c.wr_en,
      rd_en     => re.outfifo.c.valid,
      dout      => re.outfifo.dout,
      full      => re.outfifo.c.full,
      wr_ack    => re.outfifo.c.wr_ack,
      overflow  => re.outfifo.c.overflow,
      empty     => re.outfifo.c.empty,
      valid     => re.outfifo.c.valid,
      underflow => re.outfifo.c.underflow
    );

outfifo_2: entity xil_defaultlib.kernel_to_streaming_fifo
    port map (
      wr_clk    => re.clk_kernel,
      rd_clk    => i.cr.clk,
      din       => re2.outfifo.din,
      wr_en     => re2.outfifo.c.wr_en,
      rd_en     => re2.outfifo.c.valid,
      dout      => re2.outfifo.dout,
      full      => re2.outfifo.c.full,
      wr_ack    => re2.outfifo.c.wr_ack,
      overflow  => re2.outfifo.c.overflow,
      empty     => re2.outfifo.c.empty,
      valid     => re2.outfifo.c.valid,
      underflow => re2.outfifo.c.underflow
    );

  -- This fifo can potentially cross clock domains to the streaming framework.
  outfifo_3: entity xil_defaultlib.kernel_to_streaming_fifo
    port map (
      wr_clk    => re.clk_kernel,
      rd_clk    => i.cr.clk,
      din       => re3.outfifo.din,
      wr_en     => re3.outfifo.c.wr_en,
      rd_en     => re3.outfifo.c.valid,
      dout      => re3.outfifo.dout,
      full      => re3.outfifo.c.full,
      wr_ack    => re3.outfifo.c.wr_ack,
      overflow  => re3.outfifo.c.overflow,
      empty     => re3.outfifo.c.empty,
      valid     => re3.outfifo.c.valid,
      underflow => re3.outfifo.c.underflow
    );

outfifo_4: entity xil_defaultlib.kernel_to_streaming_fifo
    port map (
      wr_clk    => re.clk_kernel,
      rd_clk    => i.cr.clk,
      din       => re4.outfifo.din,
      wr_en     => re4.outfifo.c.wr_en,
      rd_en     => re4.outfifo.c.valid,
      dout      => re4.outfifo.dout,
      full      => re4.outfifo.c.full,
      wr_ack    => re4.outfifo.c.wr_ack,
      overflow  => re4.outfifo.c.overflow,
      empty     => re4.outfifo.c.empty,
      valid     => re4.outfifo.c.valid,
      underflow => re4.outfifo.c.underflow
    );

  -- Swap the first floats to go into the fifo first
  resultswap: for K in 0 to DMA_DATA_WIDTH / CU_RESULT_SIZE - 1 generate
    re.outdata((K+1)*CU_RESULT_SIZE-1 downto K*CU_RESULT_SIZE) <=  re.outfifo.dout(
                                                      DMA_DATA_WIDTH -  K    * CU_RESULT_SIZE - 1
                                                      downto
                                                      DMA_DATA_WIDTH - (K+1) * CU_RESULT_SIZE
                                                    );

  end generate;
  resultswap_2: for K in 0 to DMA_DATA_WIDTH / CU_RESULT_SIZE - 1 generate
    re2.outdata((K+1)*CU_RESULT_SIZE-1 downto K*CU_RESULT_SIZE) <=  re2.outfifo.dout(
                                                      DMA_DATA_WIDTH -  K    * CU_RESULT_SIZE - 1
                                                      downto
                                                      DMA_DATA_WIDTH - (K+1) * CU_RESULT_SIZE
                                                    );

  end generate;

  -- Swap the first floats to go into the fifo first
  resultswap_3: for K in 0 to DMA_DATA_WIDTH / CU_RESULT_SIZE - 1 generate
    re3.outdata((K+1)*CU_RESULT_SIZE-1 downto K*CU_RESULT_SIZE) <=  re3.outfifo.dout(
                                                      DMA_DATA_WIDTH -  K    * CU_RESULT_SIZE - 1
                                                      downto
                                                      DMA_DATA_WIDTH - (K+1) * CU_RESULT_SIZE
                                                    );

  end generate;
  resultswap_4: for K in 0 to DMA_DATA_WIDTH / CU_RESULT_SIZE - 1 generate
    re4.outdata((K+1)*CU_RESULT_SIZE-1 downto K*CU_RESULT_SIZE) <=  re4.outfifo.dout(
                                                      DMA_DATA_WIDTH -  K    * CU_RESULT_SIZE - 1
                                                      downto
                                                      DMA_DATA_WIDTH - (K+1) * CU_RESULT_SIZE
                                                    );

  end generate;

--------------------------------------------------------------------------------------------------- Feedback FIFO

  -- Output data that is written back to the memory, goes into the fifo first
  re.fbfifo.din  <= (        31 downto   0 => re.pairhmm.o.last.mids.ml,
                             63 downto  32 => re.pairhmm.o.last.mids.il,
                             95 downto  64 => re.pairhmm.o.last.mids.dl,
                            127 downto  96 => re.pairhmm.o.last.emis.distm_simi,
                            159 downto 128 => re.pairhmm.o.last.emis.distm_diff,
                            191 downto 160 => re.pairhmm.o.last.tmis.alpha,
                            223 downto 192 => re.pairhmm.o.last.tmis.beta,
                            255 downto 224 => re.pairhmm.o.last.tmis.delta,
                            287 downto 256 => re.pairhmm.o.last.tmis.epsilon,
                            319 downto 288 => re.pairhmm.o.last.tmis.zeta,
                            351 downto 320 => re.pairhmm.o.last.tmis.eta,
                            354 downto 352 => bpslv3(re.pairhmm.o.last.x),
                            386 downto 355 => re.pairhmm.o.last.initial,
                                    others => '0'
                          );
-- Output data that is written back to the memory, goes into the fifo first
  re2.fbfifo.din  <= (        31 downto   0 => re2.pairhmm.o.last.mids.ml,
                             63 downto  32 => re2.pairhmm.o.last.mids.il,
                             95 downto  64 => re2.pairhmm.o.last.mids.dl,
                            127 downto  96 => re2.pairhmm.o.last.emis.distm_simi,
                            159 downto 128 => re2.pairhmm.o.last.emis.distm_diff,
                            191 downto 160 => re2.pairhmm.o.last.tmis.alpha,
                            223 downto 192 => re2.pairhmm.o.last.tmis.beta,
                            255 downto 224 => re2.pairhmm.o.last.tmis.delta,
                            287 downto 256 => re2.pairhmm.o.last.tmis.epsilon,
                            319 downto 288 => re2.pairhmm.o.last.tmis.zeta,
                            351 downto 320 => re2.pairhmm.o.last.tmis.eta,
                            354 downto 352 => bpslv3(re2.pairhmm.o.last.x),
                            386 downto 355 => re2.pairhmm.o.last.initial,
                                    others => '0'
                          );


  -- Output data that is written back to the memory, goes into the fifo first
  re3.fbfifo.din  <= (        31 downto   0 => re3.pairhmm.o.last.mids.ml,
                             63 downto  32 => re3.pairhmm.o.last.mids.il,
                             95 downto  64 => re3.pairhmm.o.last.mids.dl,
                            127 downto  96 => re3.pairhmm.o.last.emis.distm_simi,
                            159 downto 128 => re3.pairhmm.o.last.emis.distm_diff,
                            191 downto 160 => re3.pairhmm.o.last.tmis.alpha,
                            223 downto 192 => re3.pairhmm.o.last.tmis.beta,
                            255 downto 224 => re3.pairhmm.o.last.tmis.delta,
                            287 downto 256 => re3.pairhmm.o.last.tmis.epsilon,
                            319 downto 288 => re3.pairhmm.o.last.tmis.zeta,
                            351 downto 320 => re3.pairhmm.o.last.tmis.eta,
                            354 downto 352 => bpslv3(re3.pairhmm.o.last.x),
                            386 downto 355 => re3.pairhmm.o.last.initial,
                                    others => '0'
                          );
-- Output data that is written back to the memory, goes into the fifo first
  re4.fbfifo.din  <= (        31 downto   0 => re4.pairhmm.o.last.mids.ml,
                             63 downto  32 => re4.pairhmm.o.last.mids.il,
                             95 downto  64 => re4.pairhmm.o.last.mids.dl,
                            127 downto  96 => re4.pairhmm.o.last.emis.distm_simi,
                            159 downto 128 => re4.pairhmm.o.last.emis.distm_diff,
                            191 downto 160 => re4.pairhmm.o.last.tmis.alpha,
                            223 downto 192 => re4.pairhmm.o.last.tmis.beta,
                            255 downto 224 => re4.pairhmm.o.last.tmis.delta,
                            287 downto 256 => re4.pairhmm.o.last.tmis.epsilon,
                            319 downto 288 => re4.pairhmm.o.last.tmis.zeta,
                            351 downto 320 => re4.pairhmm.o.last.tmis.eta,
                            354 downto 352 => bpslv3(re4.pairhmm.o.last.x),
                            386 downto 355 => re4.pairhmm.o.last.initial,
                                    others => '0'
                          );

  -- latency of 1 to match delay of read and hapl rams
  re.fbfifo.c.rd_en           <= rs.feedback_rd_en;  
  re2.fbfifo.c.rd_en          <= rs2.feedback_rd_en;
  re3.fbfifo.c.rd_en          <= rs3.feedback_rd_en;  
  re4.fbfifo.c.rd_en          <= rs4.feedback_rd_en;

  re.fbfifo.c.wr_en           <= rs.feedback_wr_en and re.pairhmm.o.last.valid;
  re2.fbfifo.c.wr_en          <= rs2.feedback_wr_en and re2.pairhmm.o.last.valid;
  re3.fbfifo.c.wr_en          <= rs3.feedback_wr_en and re3.pairhmm.o.last.valid;
  re4.fbfifo.c.wr_en          <= rs4.feedback_wr_en and re4.pairhmm.o.last.valid;

  fbfifo: entity xil_defaultlib.feedback_fifo port map (
      din       => re.fbfifo.din,
      dout      => re.fbfifo.dout,
      clk       => re.clk_kernel,
      srst      => re.fbfifo.c.rst,
      wr_en     => re.fbfifo.c.wr_en,
      rd_en     => re.fbfifo.c.rd_en,
      wr_ack    => re.fbfifo.c.wr_ack,
      valid     => re.fbfifo.c.valid,
      full      => re.fbfifo.c.full,
      empty     => re.fbfifo.c.empty,
      overflow  => re.fbfifo.c.overflow,
      underflow => re.fbfifo.c.underflow
    );

  fbfifo_2: entity xil_defaultlib.feedback_fifo port map (
      din       => re2.fbfifo.din,
      dout      => re2.fbfifo.dout,
      clk       => re.clk_kernel,
      srst      => re2.fbfifo.c.rst,
      wr_en     => re2.fbfifo.c.wr_en,
      rd_en     => re2.fbfifo.c.rd_en,
      wr_ack    => re2.fbfifo.c.wr_ack,
      valid     => re2.fbfifo.c.valid,
      full      => re2.fbfifo.c.full,
      empty     => re2.fbfifo.c.empty,
      overflow  => re2.fbfifo.c.overflow,
      underflow => re2.fbfifo.c.underflow
    );

  fbfifo_3: entity xil_defaultlib.feedback_fifo port map (
      din       => re3.fbfifo.din,
      dout      => re3.fbfifo.dout,
      clk       => re.clk_kernel,
      srst      => re3.fbfifo.c.rst,
      wr_en     => re3.fbfifo.c.wr_en,
      rd_en     => re3.fbfifo.c.rd_en,
      wr_ack    => re3.fbfifo.c.wr_ack,
      valid     => re3.fbfifo.c.valid,
      full      => re3.fbfifo.c.full,
      empty     => re3.fbfifo.c.empty,
      overflow  => re3.fbfifo.c.overflow,
      underflow => re3.fbfifo.c.underflow
    );

  fbfifo_4: entity xil_defaultlib.feedback_fifo port map (
      din       => re4.fbfifo.din,
      dout      => re4.fbfifo.dout,
      clk       => re.clk_kernel,
      srst      => re4.fbfifo.c.rst,
      wr_en     => re4.fbfifo.c.wr_en,
      rd_en     => re4.fbfifo.c.rd_en,
      wr_ack    => re4.fbfifo.c.wr_ack,
      valid     => re4.fbfifo.c.valid,
      full      => re4.fbfifo.c.full,
      empty     => re4.fbfifo.c.empty,
      overflow  => re4.fbfifo.c.overflow,
      underflow => re4.fbfifo.c.underflow
    );
  re.fbfifo.c.rst  <= rs.feedback_rst;
  re2.fbfifo.c.rst <= rs2.feedback_rst;
  re3.fbfifo.c.rst <= rs3.feedback_rst;
  re4.fbfifo.c.rst <= rs4.feedback_rst;

-- Set top left input to 1.0 when this is the first cycle of this pair.
  -- TODO: change this to be re.pairhmm_ins(J).first.mids.dtl <= a value loaded through MMIO?
  with rs.cycle select re.fbpairhmm.mids.mtl <= X"3f800000" when CYCLE_ZERO,
                                                X"00000000" when others;

  re.fbpairhmm.mids.itl         <= (others => '0');
  re.fbpairhmm.mids.dtl         <= (others => '0');
  re.fbpairhmm.mids.mt          <= (others => '0');
  re.fbpairhmm.mids.it          <= (others => '0');
  re.fbpairhmm.mids.dt          <= (others => '0');
    
  re.fbpairhmm.mids.ml          <= re.fbfifo.dout( 31 downto   0);
  re.fbpairhmm.mids.il          <= re.fbfifo.dout( 63 downto  32);
  re.fbpairhmm.mids.dl          <= re.fbfifo.dout( 95 downto  64);
  re.fbpairhmm.emis.distm_simi  <= re.fbfifo.dout(127 downto  96);
  re.fbpairhmm.emis.distm_diff  <= re.fbfifo.dout(159 downto 128);
  re.fbpairhmm.tmis.alpha       <= re.fbfifo.dout(191 downto 160);
  re.fbpairhmm.tmis.beta        <= re.fbfifo.dout(223 downto 192);
  re.fbpairhmm.tmis.delta       <= re.fbfifo.dout(255 downto 224);
  re.fbpairhmm.tmis.epsilon     <= re.fbfifo.dout(287 downto 256);
  re.fbpairhmm.tmis.zeta        <= re.fbfifo.dout(319 downto 288);
  re.fbpairhmm.tmis.eta         <= re.fbfifo.dout(351 downto 320);
  re.fbpairhmm.x                <= slv3bp(re.fbfifo.dout(354 downto 352));
  re.fbpairhmm.initial          <= re.fbfifo.dout(386 downto 355);
    
  
  re.fbpairhmm.en               <= '1';
  re.fbpairhmm.valid            <= rs.valid;
  re.fbpairhmm.cell             <= rs.cell;

-------------------------------_SECOND_------------------------------------
with rs2.cycle select re2.fbpairhmm.mids.mtl <= X"3f800000" when CYCLE_ZERO,
                                                X"00000000" when others;

  re2.fbpairhmm.mids.itl         <= (others => '0');
  re2.fbpairhmm.mids.dtl         <= (others => '0');
  re2.fbpairhmm.mids.mt          <= (others => '0');
  re2.fbpairhmm.mids.it          <= (others => '0');
  re2.fbpairhmm.mids.dt          <= (others => '0');
    
  re2.fbpairhmm.mids.ml          <= re2.fbfifo.dout( 31 downto   0);
  re2.fbpairhmm.mids.il          <= re2.fbfifo.dout( 63 downto  32);
  re2.fbpairhmm.mids.dl          <= re2.fbfifo.dout( 95 downto  64);
  re2.fbpairhmm.emis.distm_simi  <= re2.fbfifo.dout(127 downto  96);
  re2.fbpairhmm.emis.distm_diff  <= re2.fbfifo.dout(159 downto 128);
  re2.fbpairhmm.tmis.alpha       <= re2.fbfifo.dout(191 downto 160);
  re2.fbpairhmm.tmis.beta        <= re2.fbfifo.dout(223 downto 192);
  re2.fbpairhmm.tmis.delta       <= re2.fbfifo.dout(255 downto 224);
  re2.fbpairhmm.tmis.epsilon     <= re2.fbfifo.dout(287 downto 256);
  re2.fbpairhmm.tmis.zeta        <= re2.fbfifo.dout(319 downto 288);
  re2.fbpairhmm.tmis.eta         <= re2.fbfifo.dout(351 downto 320);
  re2.fbpairhmm.x                <= slv3bp(re2.fbfifo.dout(354 downto 352));
  re2.fbpairhmm.initial          <= re2.fbfifo.dout(386 downto 355);
   

  re2.fbpairhmm.en               <= '1';
  re2.fbpairhmm.valid            <= rs2.valid;
  re2.fbpairhmm.cell             <= rs2.cell;
-------------------------------THIRD------------------------------------
with rs3.cycle select re3.fbpairhmm.mids.mtl <= X"3f800000" when CYCLE_ZERO,
                                                X"00000000" when others;

  re3.fbpairhmm.mids.itl         <= (others => '0');
  re3.fbpairhmm.mids.dtl         <= (others => '0');
  re3.fbpairhmm.mids.mt          <= (others => '0');
  re3.fbpairhmm.mids.it          <= (others => '0');
  re3.fbpairhmm.mids.dt          <= (others => '0');
    
  re3.fbpairhmm.mids.ml          <= re3.fbfifo.dout( 31 downto   0);
  re3.fbpairhmm.mids.il          <= re3.fbfifo.dout( 63 downto  32);
  re3.fbpairhmm.mids.dl          <= re3.fbfifo.dout( 95 downto  64);
  re3.fbpairhmm.emis.distm_simi  <= re3.fbfifo.dout(127 downto  96);
  re3.fbpairhmm.emis.distm_diff  <= re3.fbfifo.dout(159 downto 128);
  re3.fbpairhmm.tmis.alpha       <= re3.fbfifo.dout(191 downto 160);
  re3.fbpairhmm.tmis.beta        <= re3.fbfifo.dout(223 downto 192);
  re3.fbpairhmm.tmis.delta       <= re3.fbfifo.dout(255 downto 224);
  re3.fbpairhmm.tmis.epsilon     <= re3.fbfifo.dout(287 downto 256);
  re3.fbpairhmm.tmis.zeta        <= re3.fbfifo.dout(319 downto 288);
  re3.fbpairhmm.tmis.eta         <= re3.fbfifo.dout(351 downto 320);
  re3.fbpairhmm.x                <= slv3bp(re3.fbfifo.dout(354 downto 352));
  re3.fbpairhmm.initial          <= re3.fbfifo.dout(386 downto 355);
   

  re3.fbpairhmm.en               <= '1';
  re3.fbpairhmm.valid            <= rs3.valid;
  re3.fbpairhmm.cell             <= rs3.cell;
-------------------------------FOURTH------------------------------------
with rs4.cycle select re4.fbpairhmm.mids.mtl <= X"3f800000" when CYCLE_ZERO,
                                                X"00000000" when others;

  re4.fbpairhmm.mids.itl         <= (others => '0');
  re4.fbpairhmm.mids.dtl         <= (others => '0');
  re4.fbpairhmm.mids.mt          <= (others => '0');
  re4.fbpairhmm.mids.it          <= (others => '0');
  re4.fbpairhmm.mids.dt          <= (others => '0');
    
  re4.fbpairhmm.mids.ml          <= re4.fbfifo.dout( 31 downto   0);
  re4.fbpairhmm.mids.il          <= re4.fbfifo.dout( 63 downto  32);
  re4.fbpairhmm.mids.dl          <= re4.fbfifo.dout( 95 downto  64);
  re4.fbpairhmm.emis.distm_simi  <= re4.fbfifo.dout(127 downto  96);
  re4.fbpairhmm.emis.distm_diff  <= re4.fbfifo.dout(159 downto 128);
  re4.fbpairhmm.tmis.alpha       <= re4.fbfifo.dout(191 downto 160);
  re4.fbpairhmm.tmis.beta        <= re4.fbfifo.dout(223 downto 192);
  re4.fbpairhmm.tmis.delta       <= re4.fbfifo.dout(255 downto 224);
  re4.fbpairhmm.tmis.epsilon     <= re4.fbfifo.dout(287 downto 256);
  re4.fbpairhmm.tmis.zeta        <= re4.fbfifo.dout(319 downto 288);
  re4.fbpairhmm.tmis.eta         <= re4.fbfifo.dout(351 downto 320);
  re4.fbpairhmm.x                <= slv3bp(re4.fbfifo.dout(354 downto 352));
  re4.fbpairhmm.initial          <= re4.fbfifo.dout(386 downto 355);
   

  re4.fbpairhmm.en               <= '1';
  re4.fbpairhmm.valid            <= rs4.valid;
  re4.fbpairhmm.cell             <= rs4.cell;

---------------------------------------------------------------------------------------------------
--     _____           _        _ _
--    / ____|         | |      | (_)          /\
--   | (___  _   _ ___| |_ ___ | |_  ___     /  \   _ __ _ __ __ _ _   _
--    \___ \| | | / __| __/ _ \| | |/ __|   / /\ \ | '__| '__/ _` | | | |
--    ____) | |_| \__ \ || (_) | | | (__   / ____ \| |  | | | (_| | |_| |
--   |_____/ \__, |___/\__\___/|_|_|\___| /_/    \_\_|  |_|  \__,_|\__, |
--            __/ |                                                 __/ |
--           |___/                                                 |___/
---------------------------------------------------------------------------------------------------
  -- Connect clock and reset
  re.pairhmm_cr                 <=  ( clk => re.clk_kernel,
                                      rst => rs.pairhmm_rst
                                    );

  re2.pairhmm_cr                 <=  ( clk => re.clk_kernel,
                                      rst => rs2.pairhmm_rst
                                    );

  re3.pairhmm_cr                 <=  ( clk => re.clk_kernel,
                                      rst => rs3.pairhmm_rst
                                    );

  re4.pairhmm_cr                 <=  ( clk => re.clk_kernel,
                                      rst => rs4.pairhmm_rst
                                    );

  -- Input for the first PE
  re.pairhmm.i.first            <= rs.pe_first;
  re2.pairhmm.i.first           <= rs2.pe_first;
  re3.pairhmm.i.first           <= rs3.pe_first;
  re4.pairhmm.i.first           <= rs4.pe_first;

  -- Base X for the first PE must come from the read RAM or it must come from the feedback FIFO with a latency of 1
  re.pairhmm.i.x                <= slv3bp(re.readram.doutb((idx(rs.core_schedule) + 1)* BP_SIZE -1 downto idx(rs.core_schedule)*BP_SIZE)) when rs.feedback_rd_en1 = '0' else
                                   rs.pe_first.x;

  re2.pairhmm.i.x                <= slv3bp(re2.readram.doutb((idx(rs2.core_schedule) + 1)* BP_SIZE -1 downto idx(rs2.core_schedule)*BP_SIZE)) when rs2.feedback_rd_en1 = '0' else
                                   rs2.pe_first.x;

  re3.pairhmm.i.x                <= slv3bp(re3.readram.doutb((idx(rs3.core_schedule) + 1)* BP_SIZE -1 downto idx(rs3.core_schedule)*BP_SIZE)) when rs3.feedback_rd_en1 = '0' else
                                   rs3.pe_first.x;

  re4.pairhmm.i.x                <= slv3bp(re4.readram.doutb((idx(rs4.core_schedule) + 1)* BP_SIZE -1 downto idx(rs4.core_schedule)*BP_SIZE)) when rs4.feedback_rd_en1 = '0' else
                                   rs4.pe_first.x;
  -- Schedule
  re.pairhmm.i.schedule         <= rs.core_schedule;
  re2.pairhmm.i.schedule        <= rs2.core_schedule;
  re3.pairhmm.i.schedule        <= rs3.core_schedule;
  re4.pairhmm.i.schedule        <= rs4.core_schedule;

  -- Address for Y bus
  re.pairhmm.i.ybus.addr         <= rs.ybus_addr1;
  re.pairhmm.i.ybus.wren         <= rs.ybus_en1;

  re2.pairhmm.i.ybus.addr        <= rs2.ybus_addr1;
  re2.pairhmm.i.ybus.wren        <= rs2.ybus_en1;

  re3.pairhmm.i.ybus.addr        <= rs3.ybus_addr1;
  re3.pairhmm.i.ybus.wren        <= rs3.ybus_en1;

  re4.pairhmm.i.ybus.addr        <= rs4.ybus_addr1;
  re4.pairhmm.i.ybus.wren        <= rs4.ybus_en1;
  
  re.pairhmm.i.fb                <= rs.feedback_wr_en;
  re2.pairhmm.i.fb               <= rs2.feedback_wr_en;
  re3.pairhmm.i.fb               <= rs3.feedback_wr_en;
  re4.pairhmm.i.fb               <= rs4.feedback_wr_en;
  -- Data for Y bus
  ybus_data_sel: for J in 0 to PE_DEPTH - 1 generate
    re.pairhmm.i.ybus.data(J)   <= slv3bp(re.haplram.doutb((J+1) * BP_SIZE - 1 downto J * BP_SIZE)) when rs.ybus_addr1 < rs.sizey1 else
                                   BP_STOP;
  end generate;

 ybus_data_sel_2: for J in 0 to PE_DEPTH - 1 generate
    re2.pairhmm.i.ybus.data(J)   <= slv3bp(re2.haplram.doutb((J+1) * BP_SIZE - 1 downto J * BP_SIZE)) when rs2.ybus_addr1 < rs2.sizey1 else
                                   BP_STOP;
  end generate;
  -- Data for Y bus
  ybus_data_sel_3: for J in 0 to PE_DEPTH - 1 generate
    re3.pairhmm.i.ybus.data(J)   <= slv3bp(re3.haplram.doutb((J+1) * BP_SIZE - 1 downto J * BP_SIZE)) when rs3.ybus_addr1 < rs3.sizey1 else
                                   BP_STOP;
  end generate;

 ybus_data_sel_4: for J in 0 to PE_DEPTH - 1 generate
    re4.pairhmm.i.ybus.data(J)   <= slv3bp(re4.haplram.doutb((J+1) * BP_SIZE - 1 downto J * BP_SIZE)) when rs4.ybus_addr1 < rs4.sizey1 else
                                   BP_STOP;
  end generate;
  -- Core instantiation
  pairhmm_core : entity work.pairhmm port map (
    cr  => re.pairhmm_cr,
    i   => re.pairhmm.i,
    o   => re.pairhmm.o
  );

  pairhmm_core_2 : entity work.pairhmm port map (
    cr  => re2.pairhmm_cr,
    i   => re2.pairhmm.i,
    o   => re2.pairhmm.o
  );

  pairhmm_core_3 : entity work.pairhmm port map (
    cr  => re3.pairhmm_cr,
    i   => re3.pairhmm.i,
    o   => re3.pairhmm.o
  );

  pairhmm_core_4 : entity work.pairhmm port map (
    cr  => re4.pairhmm_cr,
    i   => re4.pairhmm.i,
    o   => re4.pairhmm.o
  );
---------------------------------------------------------------------------------------------------
--     _____      _              _       _
--    / ____|    | |            | |     | |
--   | (___   ___| |__   ___  __| |_   _| | ___ _ __
--    \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \ '__|
--    ____) | (__| | | |  __/ (_| | |_| | |  __/ |
--   |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|_|
---------------------------------------------------------------------------------------------------
-- This implements a round-robin scheduler to fill pipeline stage n with the output of FIFO n
---------------------------------------------------------------------------------------------------

  scheduler_comb : process(all)
    variable vs,vs2,vs3,vs4    : cu_sched;
  begin
--------------------------------------------------------------------------------------------------- default assignments
    vs                                  := rs;
    vs2					                := rs2;
    vs3                                 := rs3;
    vs4					                := rs4;
    vs.ybus_en                          := '0';
    vs2.ybus_en                         := '0';
    vs3.ybus_en                         := '0';
    vs4.ybus_en                         := '0';
    -- Select the proper input, also correct for latency of 1 of the hapl and read RAMs:
    if rs.feedback_rd_en = '0' then
      vs.pe_first                       := re.pairhmm_in;
    else
      vs.pe_first                       := re.fbpairhmm;
    end if;

    if rs2.feedback_rd_en = '0' then
      vs2.pe_first                       := re2.pairhmm_in;
    else
      vs2.pe_first                       := re2.fbpairhmm;
    end if;

    if rs3.feedback_rd_en = '0' then
      vs3.pe_first                       := re3.pairhmm_in;
    else
      vs3.pe_first                       := re3.fbpairhmm;
    end if;

    if rs4.feedback_rd_en = '0' then
      vs4.pe_first                       := re4.pairhmm_in;
    else
      vs4.pe_first                       := re4.fbpairhmm;
    end if;

    -- Control signals that also need a latency of 1 for this reason:
    vs.ybus_addr1                        := rs.ybus_addr;
    vs.core_schedule                     := rs.schedule;
    vs.feedback_rd_en1                   := rs.feedback_rd_en;
    vs.ybus_en1                          := rs.ybus_en;
    vs.sizey1                            := rs.sizey;

    vs2.ybus_addr1                       := rs2.ybus_addr;
    vs2.core_schedule                    := rs2.schedule;
    vs2.feedback_rd_en1                  := rs2.feedback_rd_en;
    vs2.ybus_en1                         := rs2.ybus_en;
    vs2.sizey1                           := rs2.sizey;

    vs3.ybus_addr1                       := rs3.ybus_addr;
    vs3.core_schedule                    := rs3.schedule;
    vs3.feedback_rd_en1                  := rs3.feedback_rd_en;
    vs3.ybus_en1                         := rs3.ybus_en;
    vs3.sizey1                           := rs3.sizey;

    vs4.ybus_addr1                       := rs4.ybus_addr;
    vs4.core_schedule                    := rs4.schedule;
    vs4.feedback_rd_en1                  := rs4.feedback_rd_en;
    vs4.ybus_en1                         := rs4.ybus_en;
    vs4.sizey1                           := rs4.sizey;

--------------------------------------------------------------------------------------------------- round robin schedule
    -- Schedule is always running, this is to keep the PairHMM core running even when the scheduler
    -- itself is idle. This allows the scheduler to start a new batch while there is still an old
    -- batch somewhere in the Systolic Array

    -- Go to next pair
    vs.schedule                         := rs.schedule + 1;
	vs2.schedule                        := rs2.schedule + 1;
    vs3.schedule                        := rs3.schedule + 1;
	vs4.schedule                        := rs4.schedule + 1;
    -- Wrap around (required when log2(PE_DEPTH) is not an integer
    
    if vs.schedule = PE_DEPTH then
      vs.schedule                       := (others => '0');
    end if;

    if vs2.schedule = PE_DEPTH then
      vs2.schedule                       := (others => '0');
    end if;

    if vs3.schedule = PE_DEPTH then
      vs3.schedule                       := (others => '0');
    end if;

    if vs4.schedule = PE_DEPTH then
      vs4.schedule                       := (others => '0');
    end if;
--------------------------------------------------------------------------------------------------- state machine
case rs.state is

      when idle =>
        -- Gather the sizes, bases and initial D row value from the other clock domain
        vs.sizey                        := r.sched.y_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs.sizex                        := r.sched.x_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        vs.sizexp                       := r.sched.x_padded(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs.initial_array                := r.sched_array;

        vs.cycle                        := (others => '0');
        vs.basepair                     := (others => '0');
        vs.element                      := (others => '0');
        vs.supercolumn                  := (others => '0');

        vs.valid                        := '0';
        vs.cell                         := PE_NORMAL;

        vs.fifo_rd_en                   := '0';
        vs.fifo_reads                   := (others => '0');

        vs.ybus_addr                    := (others => '0');

        vs.feedback_rd_en               := '0';
        vs.feedback_wr_en               := '0';
        vs.feedback_rst                 := '1';

        -- Start everything when the FIFO's are filled and align with the scheduler
        -- Starting up takes two cycles, thus we wait until the scheduler is at PE_DEPTH - 2.
        if r.filled = '1' and rs.schedule = PE_DEPTH - 2 then
          vs.state                      := startup;
          vs.feedback_rst               := '0';
          vs.pairhmm_rst                := '0';
        end if;

      when startup =>
        vs.state                        := processing;                            -- Go to processing state
        vs.fifo_rd_en                   := '1';                                   -- Put the read-enable high for the first FIFO
        vs.valid                        := '1';                                   -- Enable data valid on the next cycle
        vs.cell                         := PE_TOP;                                -- First cycle we will be at top of matrix
        vs.startflag                    := '1';
        vs.ram                          := not rs.ram;                            -- Wrap to the other half of the RAM
        vs.ybus_en                      := '1';

      when processing =>

        -- Unset the startflag
        vs.startflag                    := '0';

        -- Increase PairHMM cycle, except when this is the first cycle, which we can check by looking at the last bit of fifo_rd_en
        -- Everything inside this if statement is triggered when a new cell update cycle starts
        if vs.schedule = u(0,PE_DEPTH_BITS) and rs.startflag = '0' then
          vs.cycle                      := rs.cycle + 1;                          -- Increase the total cycle counter
          vs.basepair                   := rs.basepair + 1;                       -- Increase the counter of Y

          if rs.element /= PAIRHMM_NUM_PES - 1 then
            vs.element                  := rs.element + 1;                        -- Increase processing highest active element in supercolumn counter
            vs.ybus_addr                := rs.ybus_addr + 1;                      -- Increase X Bus address
            vs.ybus_en                  := '1';                                   -- Write to next element
          end if;

          -- If we are done with the last padded base of X
          if vs.basepair = rs.sizexp then
            --If this is not the last base
            if rs.cell /= PE_LAST then
              vs.supercolumn            := rs.supercolumn + 1;                    -- Advance to the next supercolumn
              vs.element                := (others => '0');                       -- Reset the highest active element in supercolumn counter
              vs.ybus_addr              := (others => '0');
              vs.basepair               := (others => '0');                       -- Reset the basepair counter of Y
              vs.sizey                  := rs.sizey - PAIRHMM_NUM_PES;            -- Subtract size in the X direction
              vs.ybus_en                := '1';                                   -- Write to first element in next cycle
            end if;
            vs.fifo_rd_en               := '0';                                   -- Stop reading from the FIFO
          end if;

          -- Default PE cell state is normal:
          vs.cell                       := PE_NORMAL;

          -- If the vertical basepair we're working on is 0, we are at the top of the matrix
          if vs.basepair = 0 then
            vs.cell                     := PE_TOP;
          end if;

          -- If we are at the last base of the read
          if vs.basepair = rs.sizex - 1 then
            vs.cell                     := PE_BOTTOM;                             -- Assert the "bottom" signal
            if rs.sizey <= PAIRHMM_NUM_PES then
              vs.cell                   := PE_LAST;                               -- Assert the "last" signal
            end if;
          end if;

          -- If we fed the last base of the whole pair in the previous cell update cycle
          if rs.cell = PE_LAST then
            vs.valid                    := '0';                                   -- Next inputs are not valid anymore
            vs.cell                     := PE_NORMAL;                             -- Not last anymore TODO: this is a don't care?
            vs.state                    := done;
          end if;
        end if;

        -- Enable feedback FIFO writing when we passed the number of PE's the first time.
        if vs.cycle = PAIRHMM_NUM_PES then
            vs.feedback_wr_en           := '1';
        end if;
                        
        -- Enable feedback FIFO reading when we passed the number of padded bases the first time.
        if vs.cycle = rs.sizexp then
          vs.feedback_rd_en             := '1';
        end if;

        -- Keep track of howmany reads we've done from the FIFO
        if vs.fifo_rd_en = '1' then
          vs.fifo_reads                 := rs.fifo_reads + 1;
        end if;

      when done =>
        vs.state                        := idle;
        --vs.pairhmm_rst                  := '1'; -- Dont reset the PairHMM core here, because even after emptying the FIFO's, there may still be
                                                  -- pairs in the systolic array
        -- Reset the feedback FIFO
        vs.feedback_rd_en               := '0';
        vs.feedback_wr_en               := '0';
        vs.feedback_rst                 := '1';
    --    start<='1';
      when others =>
        null;

    end case;
  
    case rs2.state is

      when idle =>
        -- Gather the sizes, bases and initial D row value from the other clock domain
        vs2.sizey                        := r.sched_2.y_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs2.sizex                        := r.sched_2.x_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        vs2.sizexp                       := r.sched_2.x_padded(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs2.initial_array                := r.sched_array_2;

        vs2.cycle                        := (others => '0');
        vs2.basepair                     := (others => '0');
        vs2.element                      := (others => '0');
        vs2.supercolumn                  := (others => '0');

        vs2.valid                        := '0';
        vs2.cell                         := PE_NORMAL;

        vs2.fifo_rd_en                   := '0';
        vs2.fifo_reads                   := (others => '0');

        vs2.ybus_addr                    := (others => '0');

        vs2.feedback_rd_en               := '0';
        vs2.feedback_wr_en               := '0';
        vs2.feedback_rst                 := '1';

        -- Start everything when the FIFO's are filled and align with the scheduler
        -- Starting up takes two cycles, thus we wait until the scheduler is at PE_DEPTH - 2.
        if r.filled_2 = '1' and rs2.schedule = PE_DEPTH - 2 then
          vs2.state                      := startup;
          vs2.feedback_rst               := '0';
          vs2.pairhmm_rst                := '0';
        end if;

      when startup =>
        vs2.state                        := processing;                            -- Go to processing state
        vs2.fifo_rd_en                   := '1';                                   -- Put the read-enable high for the first FIFO
        vs2.valid                        := '1';                                   -- Enable data valid on the next cycle
        vs2.cell                         := PE_TOP;                                -- First cycle we will be at top of matrix
        vs2.startflag                    := '1';
        vs2.ram                          := not rs2.ram;                            -- Wrap to the other half of the RAM
        vs2.ybus_en                      := '1';

      when processing =>

        -- Unset the startflag
        vs2.startflag                    := '0';

        -- Increase PairHMM cycle, except when this is the first cycle, which we can check by looking at the last bit of fifo_rd_en
        -- Everything inside this if statement is triggered when a new cell update cycle starts
        if vs2.schedule = u(0,PE_DEPTH_BITS) and rs2.startflag = '0' then
          vs2.cycle                      := rs2.cycle + 1;                          -- Increase the total cycle counter
          vs2.basepair                   := rs2.basepair + 1;                       -- Increase the counter of Y

          if rs2.element /= PAIRHMM_NUM_PES - 1 then
            vs2.element                  := rs2.element + 1;                        -- Increase processing highest active element in supercolumn counter
            vs2.ybus_addr                := rs2.ybus_addr + 1;                      -- Increase X Bus address
            vs2.ybus_en                  := '1';                                   -- Write to next element
          end if;

          -- If we are done with the last padded base of X
          if vs2.basepair = rs2.sizexp then
            --If this is not the last base
            if rs2.cell /= PE_LAST then
              vs2.supercolumn            := rs2.supercolumn + 1;                    -- Advance to the next supercolumn
              vs2.element                := (others => '0');                       -- Reset the highest active element in supercolumn counter
              vs2.ybus_addr              := (others => '0');
              vs2.basepair               := (others => '0');                       -- Reset the basepair counter of Y
              vs2.sizey                  := rs2.sizey - PAIRHMM_NUM_PES;            -- Subtract size in the X direction
              vs2.ybus_en                := '1';                                   -- Write to first element in next cycle
            end if;
            vs2.fifo_rd_en               := '0';                                   -- Stop reading from the FIFO
          end if;

          -- Default PE cell state is normal:
          vs2.cell                       := PE_NORMAL;

          -- If the vertical basepair we're working on is 0, we are at the top of the matrix
          if vs2.basepair = 0 then
            vs2.cell                     := PE_TOP;
          end if;

          -- If we are at the last base of the read
          if vs2.basepair = rs2.sizex - 1 then
            vs2.cell                     := PE_BOTTOM;                             -- Assert the "bottom" signal
            if rs2.sizey <= PAIRHMM_NUM_PES then
              vs2.cell                   := PE_LAST;                               -- Assert the "last" signal
            end if;
          end if;

          -- If we fed the last base of the whole pair in the previous cell update cycle
          if rs2.cell = PE_LAST then
            vs2.valid                    := '0';                                   -- Next inputs are not valid anymore
            vs2.cell                     := PE_NORMAL;                             -- Not last anymore TODO: this is a don't care?
            vs2.state                    := done;
          end if;
        end if;

        -- Enable feedback FIFO writing when we passed the number of PE's the first time.
        if vs2.cycle = PAIRHMM_NUM_PES then
            vs2.feedback_wr_en           := '1';
        end if;
                        
        -- Enable feedback FIFO reading when we passed the number of padded bases the first time.
        if vs2.cycle = rs2.sizexp then
          vs2.feedback_rd_en             := '1';
        end if;

        -- Keep track of howmany reads we've done from the FIFO
        if vs2.fifo_rd_en = '1' then
          vs2.fifo_reads                 := rs2.fifo_reads + 1;
        end if;

      when done =>
        vs2.state                        := idle;
        --vs.pairhmm_rst                  := '1'; -- Dont reset the PairHMM core here, because even after emptying the FIFO's, there may still be
                                                  -- pairs in the systolic array
        -- Reset the feedback FIFO
        vs2.feedback_rd_en               := '0';
        vs2.feedback_wr_en               := '0';
        vs2.feedback_rst                 := '1';
      when others =>
        null;

    end case;

   case rs3.state is

      when idle =>
        -- Gather the sizes, bases and initial D row value from the other clock domain
        vs3.sizey                        := r.sched_3.y_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs3.sizex                        := r.sched_3.x_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        vs3.sizexp                       := r.sched_3.x_padded(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs3.initial_array                := r.sched_array_3;

        vs3.cycle                        := (others => '0');
        vs3.basepair                     := (others => '0');
        vs3.element                      := (others => '0');
        vs3.supercolumn                  := (others => '0');

        vs3.valid                        := '0';
        vs3.cell                         := PE_NORMAL;

        vs3.fifo_rd_en                   := '0';
        vs3.fifo_reads                   := (others => '0');

        vs3.ybus_addr                    := (others => '0');

        vs3.feedback_rd_en               := '0';
        vs3.feedback_wr_en               := '0';
        vs3.feedback_rst                 := '1';

        -- Start everything when the FIFO's are filled and align with the scheduler
        -- Starting up takes two cycles, thus we wait until the scheduler is at PE_DEPTH - 2.
        if r.filled_3 = '1' and rs3.schedule = PE_DEPTH - 2 then
          vs3.state                      := startup;
          vs3.feedback_rst               := '0';
          vs3.pairhmm_rst                := '0';
        end if;

      when startup =>
        vs3.state                        := processing;                            -- Go to processing state
        vs3.fifo_rd_en                   := '1';                                   -- Put the read-enable high for the first FIFO
        vs3.valid                        := '1';                                   -- Enable data valid on the next cycle
        vs3.cell                         := PE_TOP;                                -- First cycle we will be at top of matrix
        vs3.startflag                    := '1';
        vs3.ram                          := not rs3.ram;                            -- Wrap to the other half of the RAM
        vs3.ybus_en                      := '1';

      when processing =>

        -- Unset the startflag
        vs3.startflag                    := '0';

        -- Increase PairHMM cycle, except when this is the first cycle, which we can check by looking at the last bit of fifo_rd_en
        -- Everything inside this if statement is triggered when a new cell update cycle starts
        if vs3.schedule = u(0,PE_DEPTH_BITS) and rs3.startflag = '0' then
          vs3.cycle                      := rs3.cycle + 1;                          -- Increase the total cycle counter
          vs3.basepair                   := rs3.basepair + 1;                       -- Increase the counter of Y

          if rs3.element /= PAIRHMM_NUM_PES - 1 then
            vs3.element                  := rs3.element + 1;                        -- Increase processing highest active element in supercolumn counter
            vs3.ybus_addr                := rs3.ybus_addr + 1;                      -- Increase X Bus address
            vs3.ybus_en                  := '1';                                   -- Write to next element
          end if;

          -- If we are done with the last padded base of X
          if vs3.basepair = rs3.sizexp then
            --If this is not the last base
            if rs3.cell /= PE_LAST then
              vs3.supercolumn            := rs3.supercolumn + 1;                    -- Advance to the next supercolumn
              vs3.element                := (others => '0');                       -- Reset the highest active element in supercolumn counter
              vs3.ybus_addr              := (others => '0');
              vs3.basepair               := (others => '0');                       -- Reset the basepair counter of Y
              vs3.sizey                  := rs3.sizey - PAIRHMM_NUM_PES;            -- Subtract size in the X direction
              vs3.ybus_en                := '1';                                   -- Write to first element in next cycle
            end if;
            vs3.fifo_rd_en               := '0';                                   -- Stop reading from the FIFO
          end if;

          -- Default PE cell state is normal:
          vs3.cell                       := PE_NORMAL;

          -- If the vertical basepair we're working on is 0, we are at the top of the matrix
          if vs3.basepair = 0 then
            vs3.cell                     := PE_TOP;
          end if;

          -- If we are at the last base of the read
          if vs3.basepair = rs3.sizex - 1 then
            vs3.cell                     := PE_BOTTOM;                             -- Assert the "bottom" signal
            if rs3.sizey <= PAIRHMM_NUM_PES then
              vs3.cell                   := PE_LAST;                               -- Assert the "last" signal
            end if;
          end if;

          -- If we fed the last base of the whole pair in the previous cell update cycle
          if rs3.cell = PE_LAST then
            vs3.valid                    := '0';                                   -- Next inputs are not valid anymore
            vs3.cell                     := PE_NORMAL;                             -- Not last anymore TODO: this is a don't care?
            vs3.state                    := done;
          end if;
        end if;

        -- Enable feedback FIFO writing when we passed the number of PE's the first time.
        if vs3.cycle = PAIRHMM_NUM_PES then
            vs3.feedback_wr_en           := '1';
        end if;
                        
        -- Enable feedback FIFO reading when we passed the number of padded bases the first time.
        if vs3.cycle = rs2.sizexp then
          vs3.feedback_rd_en             := '1';
        end if;

        -- Keep track of howmany reads we've done from the FIFO
        if vs3.fifo_rd_en = '1' then
          vs3.fifo_reads                 := rs2.fifo_reads + 1;
        end if;

      when done =>
        vs3.state                        := idle;
        --vs.pairhmm_rst                  := '1'; -- Dont reset the PairHMM core here, because even after emptying the FIFO's, there may still be
                                                  -- pairs in the systolic array
        -- Reset the feedback FIFO
        vs3.feedback_rd_en               := '0';
        vs3.feedback_wr_en               := '0';
        vs3.feedback_rst                 := '1';
      when others =>
        null;

    end case;

case rs4.state is

      when idle =>
        -- Gather the sizes, bases and initial D row value from the other clock domain
        vs4.sizey                        := r.sched_4.y_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs4.sizex                        := r.sched_4.x_size(log2e(PAIRHMM_MAX_SIZE) downto 0);
        vs4.sizexp                       := r.sched_4.x_padded(log2e(PAIRHMM_MAX_SIZE) downto 0);
        
        vs4.initial_array                := r.sched_array_4;

        vs4.cycle                        := (others => '0');
        vs4.basepair                     := (others => '0');
        vs4.element                      := (others => '0');
        vs4.supercolumn                  := (others => '0');

        vs4.valid                        := '0';
        vs4.cell                         := PE_NORMAL;

        vs4.fifo_rd_en                   := '0';
        vs4.fifo_reads                   := (others => '0');

        vs4.ybus_addr                    := (others => '0');

        vs4.feedback_rd_en               := '0';
        vs4.feedback_wr_en               := '0';
        vs4.feedback_rst                 := '1';

        -- Start everything when the FIFO's are filled and align with the scheduler
        -- Starting up takes two cycles, thus we wait until the scheduler is at PE_DEPTH - 2.
        if r.filled_4 = '1' and rs4.schedule = PE_DEPTH - 2 then
          vs4.state                      := startup;
          vs4.feedback_rst               := '0';
          vs4.pairhmm_rst                := '0';
        end if;

      when startup =>
        vs4.state                        := processing;                            -- Go to processing state
        vs4.fifo_rd_en                   := '1';                                   -- Put the read-enable high for the first FIFO
        vs4.valid                        := '1';                                   -- Enable data valid on the next cycle
        vs4.cell                         := PE_TOP;                                -- First cycle we will be at top of matrix
        vs4.startflag                    := '1';
        vs4.ram                          := not rs4.ram;                            -- Wrap to the other half of the RAM
        vs4.ybus_en                      := '1';

      when processing =>

        -- Unset the startflag
        vs4.startflag                    := '0';

        -- Increase PairHMM cycle, except when this is the first cycle, which we can check by looking at the last bit of fifo_rd_en
        -- Everything inside this if statement is triggered when a new cell update cycle starts
        if vs4.schedule = u(0,PE_DEPTH_BITS) and rs4.startflag = '0' then
          vs4.cycle                      := rs4.cycle + 1;                          -- Increase the total cycle counter
          vs4.basepair                   := rs4.basepair + 1;                       -- Increase the counter of Y

          if rs4.element /= PAIRHMM_NUM_PES - 1 then
            vs4.element                  := rs4.element + 1;                        -- Increase processing highest active element in supercolumn counter
            vs4.ybus_addr                := rs4.ybus_addr + 1;                      -- Increase X Bus address
            vs4.ybus_en                  := '1';                                   -- Write to next element
          end if;

          -- If we are done with the last padded base of X
          if vs4.basepair = rs4.sizexp then
            --If this is not the last base
            if rs4.cell /= PE_LAST then
              vs4.supercolumn            := rs4.supercolumn + 1;                    -- Advance to the next supercolumn
              vs4.element                := (others => '0');                       -- Reset the highest active element in supercolumn counter
              vs4.ybus_addr              := (others => '0');
              vs4.basepair               := (others => '0');                       -- Reset the basepair counter of Y
              vs4.sizey                  := rs4.sizey - PAIRHMM_NUM_PES;            -- Subtract size in the X direction
              vs4.ybus_en                := '1';                                   -- Write to first element in next cycle
            end if;
            vs4.fifo_rd_en               := '0';                                   -- Stop reading from the FIFO
          end if;

          -- Default PE cell state is normal:
          vs4.cell                       := PE_NORMAL;

          -- If the vertical basepair we're working on is 0, we are at the top of the matrix
          if vs4.basepair = 0 then
            vs4.cell                     := PE_TOP;
          end if;

          -- If we are at the last base of the read
          if vs4.basepair = rs4.sizex - 1 then
            vs4.cell                     := PE_BOTTOM;                             -- Assert the "bottom" signal
            if rs4.sizey <= PAIRHMM_NUM_PES then
              vs4.cell                   := PE_LAST;                               -- Assert the "last" signal
            end if;
          end if;

          -- If we fed the last base of the whole pair in the previous cell update cycle
          if rs4.cell = PE_LAST then
            vs4.valid                    := '0';                                   -- Next inputs are not valid anymore
            vs4.cell                     := PE_NORMAL;                             -- Not last anymore TODO: this is a don't care?
            vs4.state                    := done;
          end if;
        end if;

        -- Enable feedback FIFO writing when we passed the number of PE's the first time.
        if vs4.cycle = PAIRHMM_NUM_PES then
            vs4.feedback_wr_en           := '1';
        end if;
                        
        -- Enable feedback FIFO reading when we passed the number of padded bases the first time.
        if vs4.cycle = rs4.sizexp then
          vs4.feedback_rd_en             := '1';
        end if;

        -- Keep track of howmany reads we've done from the FIFO
        if vs4.fifo_rd_en = '1' then
          vs4.fifo_reads                 := rs4.fifo_reads + 1;
        end if;

      when done =>
        vs4.state                        := idle;
        --vs.pairhmm_rst                  := '1'; -- Dont reset the PairHMM core here, because even after emptying the FIFO's, there may still be
                                                  -- pairs in the systolic array
        -- Reset the feedback FIFO
        vs4.feedback_rd_en               := '0';
        vs4.feedback_wr_en               := '0';
        vs4.feedback_rst                 := '1';
      when others =>
        null;

    end case;
--------------------------------------------------------------------------------------------------- outputs
    qs                         <= vs;
    qs2                        <= vs2;
    qs3                        <= vs3;
    qs4                        <= vs4;
  end process;
--------------------------------------------------------------------------------------------------- registers
  scheduler_reg : process(re.clk_kernel)
  begin
    if rising_edge(re.clk_kernel) then
      if i.cr.rst then
        rs.state              <= idle;
        rs.cycle              <= (others => '0');
        rs.basepair           <= (others => '0');
        rs.schedule           <= (others => '0');
        rs.fifo_rd_en         <= '0';
        rs.valid              <= '0';
        rs.cell               <= PE_NORMAL;
        rs.pairhmm_rst        <= '1';
        rs.feedback_rd_en     <= '0';
        rs.feedback_wr_en     <= '0';
        rs.feedback_rst       <= '1';
        rs.fifo_reads         <= (others => '0');
        rs.ram                <= '1'; -- Set to use top half of RAM since at start this will wrap to bottom
        rs.ybus_en            <= '0';
        rs.ybus_en1           <= '0'; 

		rs2.state              <= idle;
        rs2.cycle              <= (others => '0');
        rs2.basepair           <= (others => '0');
        rs2.schedule           <= (others => '0');
        rs2.fifo_rd_en         <= '0';
        rs2.valid              <= '0';
        rs2.cell               <= PE_NORMAL;
        rs2.pairhmm_rst        <= '1';
        rs2.feedback_rd_en     <= '0';
        rs2.feedback_wr_en     <= '0';
        rs2.feedback_rst       <= '1';
        rs2.fifo_reads         <= (others => '0');
        rs2.ram                <= '1'; -- Set to use top half of RAM since at start this will wrap to bottom
        rs2.ybus_en            <= '0';
        rs2.ybus_en1           <= '0';   

		rs3.state              <= idle;
        rs3.cycle              <= (others => '0');
        rs3.basepair           <= (others => '0');
        rs3.schedule           <= (others => '0');
        rs3.fifo_rd_en         <= '0';
        rs3.valid              <= '0';
        rs3.cell               <= PE_NORMAL;
        rs3.pairhmm_rst        <= '1';
        rs3.feedback_rd_en     <= '0';
        rs3.feedback_wr_en     <= '0';
        rs3.feedback_rst       <= '1';
        rs3.fifo_reads         <= (others => '0');
        rs3.ram                <= '1'; -- Set to use top half of RAM since at start this will wrap to bottom
        rs3.ybus_en            <= '0';
        rs3.ybus_en1           <= '0';

	    rs4.state              <= idle;
        rs4.cycle              <= (others => '0');
        rs4.basepair           <= (others => '0');
        rs4.schedule           <= (others => '0');
        rs4.fifo_rd_en         <= '0';
        rs4.valid              <= '0';
        rs4.cell               <= PE_NORMAL;
        rs4.pairhmm_rst        <= '1';
        rs4.feedback_rd_en     <= '0';
        rs4.feedback_wr_en     <= '0';
        rs4.feedback_rst       <= '1';
        rs4.fifo_reads         <= (others => '0');
        rs4.ram                <= '1'; -- Set to use top half of RAM since at start this will wrap to bottom
        rs4.ybus_en            <= '0';
        rs4.ybus_en1           <= '0';
      else
        rs                    <= qs;
		rs2                    <= qs2;
		rs3                    <= qs3;
		rs4                    <= qs4;
      end if;
    end if;
end process;

end architecture logic;
