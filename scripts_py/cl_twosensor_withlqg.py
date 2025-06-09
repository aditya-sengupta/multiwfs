import numpy as np

def run_cl_2sensor(input_ts, decimation, delayfr=1, \
                   fast_gain=0.2, slow_gain=0.2, integ=0.995, \
                   ar1hp_coeff=None, \
                   mnoisetimeseries=None, slowmnoisetimeseries=None, \
                   ncptimeseries=None, slowncptimeseries=None, \
                   separate_lqgs=None, combined_lqg=None, \
                   debugMA=False, debugrepeat=False):
    """ run a very simple single variable CL to look at ETF modeling quality"""

    """
    Inputs:
    input_ts: an array that is the input time series (phi) that the AO system will correct
    decimation: the difference in frame rate between fast and slow, e.g. 10
    delayfr: the delay in terms of # of frames, default 1
    fast_gain: the gain on the fast WFS when IC is used.
    slow_gain: the gain on the slow WFS when IC is used.
    integ: the integrator value on the fast and slow WFSs when IC is used
    ar1hp_coeff: the coefficient of the AR(1) HPF. When None, HPF is not used (default)

    Then there are other time series that we can inject:
    mnoisetimeseries: an array that is the time series of measurement noise on the fast WFS
    slowmnoisetimeseries: an array that is the time series of measurement noise on the slow WFS
    ncptimeseries: an array that is the time series of the NCP on the fast WFS
    slowncptimeseries: an array that is the time series of the NCP on the slow WFS

    separate_lqgs: a 2-item list that has dictionaries with the LQG filter matrices/vectors.
                This is how we have LQG-first implemented.
    combined_lqg: a dictionary with the LQG filter matrices/vectors. This mode uses both the
                fast and slow WFS measurements every time step (repeating the slow for decimation
                frames) in one larger LQG model.

    the debugging flags you won't use

    The outputs are:
    output_ts: the error immediately after the DM (w/o any NCP)
    meas_slow_ts: the measurements on the slow WFS
    meas_fast_ts: the measurements on the fast WFS
    meas_fast_afterfilter_ts: the measurements on the fast WFS after an HPF has applied
    phase_seen_by_slow_sensor: the residual phase including the slow NCP (this is X)

    """


    dec = int(decimation)
    slow_buffer = np.zeros(int(delayfr)+dec)
    fast_buffer = np.zeros(int(delayfr)+1)       


    useAR1 = False
    if ar1hp_coeff is not None:
        # !! using the AR(1) HP filter instead!
        useAR1 = True
        hp_alpha = ar1hp_coeff
        print("Using alpha = %0.3f for simple HPF" %(hp_alpha))

    steps = len(input_ts)

    if ncptimeseries is not None:
        if len(ncptimeseries) != steps:
            raise Exception('error - ncptimeseries data must be same length as input_ts')
    if slowncptimeseries is not None:
        if len(slowncptimeseries) != steps:
            raise Exception('error - slowncptimeseries data must be same length as input_ts')
    if mnoisetimeseries is not None:
        if len(mnoisetimeseries) != steps:
            raise Exception('error - mnoisetimeseries data must be same length as input_ts')
    if slowmnoisetimeseries is not None:
        if len(slowmnoisetimeseries) != steps:
            raise Exception('error - slowmnoisetimeseries data must be same length as input_ts')

    output_ts = np.zeros((steps))
    meas_slow_ts = np.zeros((steps))
    meas_fast_ts = np.zeros((steps))
    meas_fast_afterfilter_ts = np.zeros((steps))
    phase_seen_by_slow_sensor = np.zeros((steps))   

    # now determine if we are using LQG!
    useCombinedLQG = False
    if combined_lqg is not None:
        useCombinedLQG = True     

        nstate = np.shape(combined_lqg['ikca'])[0]
        last_state = np.zeros((nstate, 1))       
        before_ho_commands = np.zeros((1, 1))
        this_measurement = np.zeros((combined_lqg['nmeas'], 1))

    useSeparateLQGs = False
    if separate_lqgs is not None:
        if len(separate_lqgs) != 2:
            raise Exception('separate_lqgs inpu must be a 2-item list')
        useSeparateLQGs = True

        fast_nstate = np.shape(separate_lqgs[0]['ikca'])[0]
        fast_last_state = np.zeros((fast_nstate, 1))       
        fast_before_ho_commands = np.zeros((1, 1))
        fast_this_measurement = np.zeros((1, 1))

        slow_nstate = np.shape(separate_lqgs[1]['ikca'])[0]
        slow_last_state = np.zeros((slow_nstate, 1))       
        slow_before_ho_commands = np.zeros((1, 1))
        slow_this_measurement = np.zeros((1, 1))


    if useCombinedLQG and useSeparateLQGs:
        raise Exception('you cannot select both separate_lqgs and combined_lqg!')

    useLQG = useCombinedLQG or useSeparateLQGs

    this_dm = 0.
    fast_this_c = 0.
    fast_last_c = 0.
    slow_this_s = 0.
    slow_this_c = 0.
    slow_last_c = 0.
    ar_thisx = 0.
    ar_lastx = 0.
    ar_thisy = 0.
    ar_lasty = 0.


    for t in range(steps):
        this_e = input_ts[t] - this_dm
        output_ts[t] = this_e

        # ----------------------------
        # fast sensor
        fast_buffer[-1] = this_e

        # assuming it's only on the fast sensor
        if ncptimeseries is not None:
            fast_buffer[-1] += ncptimeseries[t]

        # get the measurement, delayed
        fast_this_s = fast_buffer[0]


        # assuming it's only on the fast sensor
        if mnoisetimeseries is not None:
            fast_this_s += mnoisetimeseries[t]

        # save for later analysis - the noisy measurements seen by the computer
        meas_fast_ts[t] = fast_this_s       

        # is using, apply the HPF
        if useAR1:
            ar_thisx = fast_this_s
            ar_thisy = hp_alpha*ar_lasty + hp_alpha*(ar_thisx - ar_lastx)
            ar_lastx = ar_thisx
            ar_lasty = ar_thisy
            fast_this_s = ar_thisy

        # save for later analysis - the noisy measurements after filtering by the computer
        meas_fast_afterfilter_ts[t] = fast_this_s

        # advance FAST buffer
        fast_buffer[0:-1] = fast_buffer[1:]



        # ----------------------------
        # get the SLOW measurement, delayed
        slow_buffer[-1] = this_e

        # assumimng it's only on the fast sensor
        if slowncptimeseries is not None:
            slow_buffer[-1] += slowncptimeseries[t]

        # get the measurement, delayed
        phase_seen_by_slow_sensor[t] = slow_buffer[0]

        # do we use the (linear, time-invariant) MA every single
        # fast time step, or do we take the average only once every
        # 1/sim_decimation timesteps?
        if debugMA:
            # this will update the moving-average on the
            # slow sensor every fast time step, WHICH IS NOT PHYSICAL
            slow_this_s = slow_buffer[0]
            # slow_this_s = np.sum(slow_buffer[0:dec])/dec
        else:
            # the actual physical slow sensor, which samples every decimation steps...
            if np.mod(t, dec) == (dec-1):
                # end of frame! grab out!
                slow_this_s = np.sum(slow_buffer[0:dec])/dec

        # assuming it's only on the slow sensor
        if slowmnoisetimeseries is not None:
            if debugMA:
                slow_this_s += slowmnoisetimeseries[t]
            else:
                if np.mod(t, dec) == (dec-1):
                    slow_this_s += slowmnoisetimeseries[t]

        # save for later analysis - the noisy measurements seen by the computer
        meas_slow_ts[t] = slow_this_s   

        # advance SLOW buffer
        slow_buffer[0:-1] = slow_buffer[1:]*1.






        # Now apply the controllers!
        if useLQG:
            if useCombinedLQG:
                # combined LQG is an LQG at every (fast) frame that
                # uses both measurements at once!

                if combined_lqg['fastindex'] is not None:
                    this_measurement[combined_lqg['fastindex'],0] = fast_this_s
                if debugrepeat or debugMA:
                    if combined_lqg['slowindex'] is not None:
                        this_measurement[combined_lqg['slowindex'],0] = slow_this_s
                else:
                    if np.mod(t, dec) == (dec-1):
                        if combined_lqg['slowindex'] is not None:
                            this_measurement[combined_lqg['slowindex'],0] = slow_this_s

                state_part1 = np.dot(combined_lqg['k'], this_measurement)
                # print(np.shape(state_part1))
                # predict the state
                state_part2 = np.dot(combined_lqg['ikca'], last_state)
                # print(np.shape(state_part2))
                # combine
                state_part3 = state_part1 + state_part2
                # print(np.shape(state_part3))
                # predict the previous commands
                state_part4 = np.dot(combined_lqg['ikcd'], before_ho_commands)
                # print(np.shape(state_part4))
                # combine
                last_state = state_part3 + state_part4
                # print(np.shape(self.before_ho_state))
                # now apply second kalman gain to get actuator ocmmands
                before_ho_commands = np.dot(combined_lqg['g'], last_state)
                # print(np.shape(self.before_ho_commands))
                this_c = before_ho_commands[0][0]
                this_dm = this_c
            else:
                # this is two separate LQGs, each operating independently
                # at their own frame rate, much like the IC mode....

                if debugrepeat or debugMA:
                    raise Exception('debugrepeat and debugMA modes not supported with separateLQGs!')
                else:
                    if np.mod(t, dec) == (dec-1):
                        slow_state_part1 = np.dot(separate_lqgs[1]['k'], slow_this_s)
                        # print(np.shape(state_part1))
                        # predict the state
                        slow_state_part2 = np.dot(separate_lqgs[1]['ikca'], slow_last_state)
                        # print(np.shape(state_part2))
                        # combine
                        slow_state_part3 = slow_state_part1 + slow_state_part2
                        # print(np.shape(state_part3))
                        # predict the previous commands
                        slow_state_part4 = np.dot(separate_lqgs[1]['ikcd'], slow_before_ho_commands)
                        # print(np.shape(state_part4))
                        # combine
                        slow_last_state = slow_state_part3 + slow_state_part4
                        # print(np.shape(self.before_ho_state))
                        # now apply second kalman gain to get actuator ocmmands
                        slow_before_ho_commands = np.dot(separate_lqgs[1]['g'], slow_last_state)
                        # print(np.shape(self.before_ho_commands))
                        slow_this_c = slow_before_ho_commands[0][0]*1

                fast_state_part1 = np.dot(separate_lqgs[0]['k'], fast_this_s)
                # print(np.shape(state_part1))
                # predict the state
                fast_state_part2 = np.dot(separate_lqgs[0]['ikca'], fast_last_state)
                # print(np.shape(state_part2))
                # combine
                fast_state_part3 = fast_state_part1 + fast_state_part2
                # print(np.shape(state_part3))
                # predict the previous commands
                fast_state_part4 = np.dot(separate_lqgs[0]['ikcd'], fast_before_ho_commands)
                # print(np.shape(state_part4))
                # combine
                fast_last_state = fast_state_part3 + fast_state_part4
                # print(np.shape(self.before_ho_state))
                # now apply second kalman gain to get actuator ocmmands
                fast_before_ho_commands = np.dot(separate_lqgs[0]['g'], fast_last_state)
                # print(np.shape(self.before_ho_commands))
                fast_this_c = fast_before_ho_commands[0][0]*1

                # combine the two control outputs
                this_dm = fast_this_c + slow_this_c

        else:       
            # apply the slow sensor control law
            if debugrepeat or debugMA:
                slow_this_c = slow_gain*slow_this_s + integ*slow_last_c
                slow_last_c = slow_this_c*1.
            else:
                # actual physical approach - everyth 10th step...
                if np.mod(t, dec) == (dec-1):
                    slow_this_c = slow_gain*slow_this_s + integ*slow_last_c
                    slow_last_c = slow_this_c*1.

            # here's where we apply the control law for the fast sensor!
            fast_this_c = fast_gain*fast_this_s + integ*fast_last_c
            fast_last_c = fast_this_c*1.

            # combine the two control outputs
            this_dm = fast_this_c + slow_this_c

    # return all the tracked signals
    return output_ts, meas_slow_ts, meas_fast_ts, meas_fast_afterfilter_ts, phase_seen_by_slow_sensor
