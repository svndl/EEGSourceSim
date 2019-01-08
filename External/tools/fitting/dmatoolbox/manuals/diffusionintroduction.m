%% An introduction to the diffusion model
%% Background
% In many settings in experimental psychology, researchers collect more
% than one type of data. In fact, as the complexity of psychological
% theories grows, the number of dependent variables in experiments is bound
% to increase, and the information contained in the interrelations between
% these dependents becomes more interesting. Advances in psychological
% statistics, psychometrics and mathematical psychology should (and do)
% reflect this change in the type of data we are prone to encounter.
% However, standard analyses for multivariate data are hitherto only
% available for some rather specific situations (e.g., when all the
% dependents can be assumed to be normally distributed, multivariate linear
% models can be applied, etc.).
%%
% A somewhat common situation where these latter conditions do not apply,
% however, is the combination of reaction time and accuracy data (a
% ubiquitous combination often referred to as two-choice response time
% data). Tasks put to participants in psychological experiments often
% combine these two variables, one of which is binary (the accuracy), while
% the other is strictly positive and strongly right-skewed. Thus, both are
% decidedly non-normally distributed and classical linear models are hence
% unsuitable.
%%
% For the analysis of two-choice response time data, several nonlinear
% stochastic models have been developed, often with substantive
% interpretations attached to them (e.g., the discrete random walk model;
% Laming, 1968; Link & Heath, 1975). A more advanced model – and the one
% that is at the heart of the present article – was first proposed in 1978
% (the Ratcliff diffusion model; Ratcliff, 1978, 1981, 1985, 1987, 1988).
% The latter model, which will be described in detail in the next section,
% has performed remarkably well in the analysis of two-choice response time
% data. It has successfully been applied to experiments in many different
% fields, such as memory (Ratcliff, 1978, 1988), letter-matching (Ratcliff,
% 1981), lexical decision (Ratcliff, Gomez, & McKoon, 2004), signal
% detection (Ratcliff & Rouder, 1998; Ratcliff, Thapar, & McKoon, 2001;
% Ratcliff, Van Zandt, & McKoon, 1999), visual search (Strayer & Kramer,
% 1994), and perceptual judgment (Ratcliff, 2002; Ratcliff & Rouder, 2000;
% Thapar, Ratcliff, & McKoon, 2003).

%% The diffusion process
% The diffusion process has been used to describe and model
% the decision component in simple two-choice tasks. In the model, it is
% assumed that an observer has a unidimensional internal representation of
% evidence, where the two choices reside at the extremes. When the observer
% is presented with a stimulus, information regarding it is accumulated
% sequentially over time and the amount of evidence changes. This feature
% is shared with other sequential sampling models such as the discrete
% random walk model (Laming, 1968; Link & Heath, 1975), the accumulator
% model (Smith & Vickers, 1988; Vickers, 1970), and the Poisson counter
% model (LaBerge, 1994; Pike, 1966, 1973; Townsend & Ashby, 1983). The
% evidence is accumulated towards one or the other response criterion until
% the total exceeds one of the criteria. These two criteria are called
% absorbing boundaries, each of which is associated with one response. The
% decision time is defined as the time from the start of the process until
% the moment one of the absorbing boundaries is reached.
%% Boundary separation
% Thus, a first parameter of the Ratcliff diffusion model is the boundary
% separation, denoted by _a_. If this boundary separation is small, the
% process is expected to end sooner but it is more prone to error since
% random variability inherent to the decision process may cause it to end
% up at the wrong boundary. When the boundaries are separated further from
% each other, both accuracy and expected reaction time will increase. The
% distance between the two absorbing boundaries therefore regulates the
% relation between speed and accuracy (the so-called speed-accuracy
% trade-off).
%% Starting point
% The information accumulation process itself has three other important
% characteristics. First, the point at which the accumulation process
% starts is called the starting point, which is denoted as _z0_ here. The
% starting point may be located closer to one of the boundaries than to the
% other but it always lies between them (0 < _z0_ < _a_). This parameter
% introduces the possibility of response bias in the decision process
% because the process is more likely to end at the boundary closer to the
% starting point. 
%% Between-trial variability in starting point
% Instead of considering the starting point fixed, we will assume it to
% vary from trial to trial according to a uniform distribution, with mean
% _z_ (0 < _z0_ < _a_) and range _sz_ (o < _sz_ < min(_z_, _a-z_). The
% choice for a uniform distribution is obvious in
% this case because the starting point has to be restricted to lie between
% the two absorbing boundaries. In our case, the end points of the uniform
% distribution can be set so that z0 can never exceed either of the
% absorbing boundaries. This would not be posssible with, for example, the
% Gaussian distribution.
%% Drift rate
% A second feature is that the information accumulation process, depending
% on the stimulus presented, can have a tendency to drift off to one of the
% two absorbing boundaries. This tendency is called the drift rate, denoted
% by _xi_. The drift rate is determined by the quality of the stimulus. A
% non-ambiguous stimulus leads to a large positive drift rate such that the
% process has a high probability of hitting the upper boundary (indicating
% a correct response) in a short time. However, for a highly ambiguous
% stimulus, the drift rate will be close to zero so that there is an almost
% 50% chance to end up at either of the two boundaries. Moreover, the
% decision process takes longer when the drift rate is close to zero
% (because then there is no outspoken tendency to drift off to either
% boundary).
%% Within-trial variability in drift rate
% A third feature of the diffusion process is that the information
% accumulation process is stochastic, implying that there is some
% uncertainty in the outcome and time course of the decision process. The
% drift rate _v_  quantifies only the mean rate of information accumulation
% but the actual drift is subject to variability and this amount of
% variability is captured by the parameter _s_, representing the volatility
% of the accumulation process. This standard deviation is a non-identified
% parameter in the model because it is involved in a trade-off relation
% with the other parameters of the model. This means that the model does
% not change (i.e., the same predictions for reaction times and choice
% probabilities follow) if all parameters are multiplied by the same
% constant. Therefore, we must fix s to an arbitrary value, which will be
% 0.1 in this case because that corresponds to a consensus value in the
% literature (e.g. Ratcliff et al., 1999; Tuerlinckx, 2004). The
% variability of the information accumulation process within a trial
% implies that the precise course of the decision trajectory is
% unpredictable. Due to the inherently stochastic nature of the model,
% different decision processes can evolve differently, even if they are
% governed by the same set of parameter values. A consequence of this is
% that both the time to decide and the final choice are random variables.
%% Between-trial variability in drift rate
% Besides the deviations of the information accumulation process from the
% drift rate within a trial, it is also quite unlikely that the invoked
% drift rate _xi_ will be the same every time a similar stimulus is presented
% to a person. There are two reasons that force us to reject that
% possibility. First, stimuli may be formally the same but they will be
% seldom exactly the same on all possible dimensions. Such variation will
% lead to different drift rates _xi_. An example is a pixel discrimination
% experiment where the ratio of black to white pixels may be the same for a
% class of stimuli, but the exact configuration of black and white pixels
% differs over these stimuli with the same ratio. Secondly, even it were
% possible to present exactly the same stimulus to a person twice, the
% perceived quality of the stimulus will differ on both occasions due to
% (internal and external) perceptual noise. For these reasons we will allow
% the drift rate to differ from trial to trial, despite the fact that the
% same stimulus may be presented to the person. The variability in the
% drift rate over trials is represented by a Gaussian distribution with
% mean _v_ and standard deviation _eta_. In sum, this means that the actual rate
% of information accumulation is subject to two sources of variability:
% First, a single trial's drift rate is a draw from a normal distribution
% and, second, the rate of accumulation may fluctuate, within a trial,
% about the mean drift rate because of the assumed variability in the
% accumulation process.
%% 
% There are many conceptual reasons for assuming variability in the
% starting point and drift rate (i.e., quantities showing variability over
% trials render a-priori a more realistic model), but there are also
% important empirical arguments for taking into account these sources of
% variability. In many experiments the error reaction times can be
% systematically shorter than their correct counterparts for some
% conditions. Such phenomena could be modeled by assuming variability in
% the starting point. Error reaction times that turn out to be
% systematically longer than correct ones can be modeled by assuming
% variability in the mean drift rate. If both sources of variability are
% present in the model, even more complicated crossover patterns in the
% data (errors faster than corrects in some conditions but the reverse for
% others) can be explained (for a more detailed discussion on these topics,
% see Ratcliff & Rouder, 1998, and Ratcliff et al., 1999). Moreover,
% without the variability in drift rate the probability of a correct
% response would approach to one, and that is a highly unrealistic
% prediction, because even in the easiest conditions of an experiment, a
% small percentage of errors (lapses) is almost always observed.
%% Nondecision time
% Another component of the model we will be using in this paper represents
% that part of the observed reaction time that did not originate from the
% diffusion process (termed the residual time or non-decision time). It is
% the time needed to perform non-decision processes such as encoding of the
% stimulus, response preparation and execution of the motor response (Luce,
% 1986). We denote the non-decision part of the observed reaction time as
% _ter_.
%% Between-trial variability in nondecision time
% Most likely, the component times of _ter_ will vary from trial to
% trial, so we assume that _ter_ also varies randomly from trial to trial.
% To capture this variability in _ter_, we have chosen to work with the
% uniform distribution because of its simple form and parameters that are
% easy to interpret. The mean of the distribution equals _Ter_ (the mean
% time devoted to the processes other than deciding) and its range is _st_,
% so that the uniform density function ranges from _Ter_ - _st_ to _Ter_ +
% _st_. The density at any point within the boundaries is constant and
% equal to 1 / _st_ .
%% References
% * LaBerge, D. (1962). A recruitment theory of simple behavior.
% _Psychometrika, 27,_ 375-396.
%
% * Laming, D. (1968). _Information theory of choice reaction time._ New
% York: Wiley. 
%
% * Link, S. & Heath, R. (1975). A sequential theory of psychological
% discrimination. _Psychometrika, 40,_ 77-105. 
%
% * Luce, D. (1986). _Response times._ New York: Oxford University Press.
%
% * Pike, A. (1966). Stochastic models of choice behaviour: Response
% probabilities and latencies of finite Markov chain systems. _British
% Journal of Mathematical and Statistical Psychology, 19,_ 15-32.
%
% * Pike, A. (1973). Response latency mechanisms for signal detection.
% _Psychological Review, 80,_ 53-68. 
%
% * Ratcliff, R. (1978). A theory of memory retrieval. _Psychological Review,
% 85,_ 59-108. 
%
% * Ratcliff, R. (1981). A theory of order relations in perceptual matching.
% _Psychological Review, 88,_ 552-572. 
%
% * Ratcliff, R. (1985). Theoretical interpretations of speed and accuracy of
% positive and negative responses. _Psychological Review, 92,_ 212-225.
%
% * Ratcliff, R. (1987). More on the speed and accuracy of positive and
% negative responses. _Psychological Review, 94,_ 277-280.  
%
% * Ratcliff, R. (1988). Continuous versus discrete information processing:
% Modeling the accumulation of partial information. _Psychological Review,_
% 95, 238-255.
%
% * Ratcliff, R. (2002). A diffusion model account of reaction time and
% accuracy in a brightness discrimination task: Fitting real data and
% failing to fit fake but plausible data. _Psychonomic Bulletin & Review,
% 9,_ 278–291.
%
% * Ratcliff, R., Gomez, P., & McKoon, G. (2004). A diffusion model account
% of the lexical-decision task. _Psychological Review, 111,_ 159-182. 
%
% * Ratcliff, R., & Rouder, J.N. (1998). Modeling response times for
% two-choice decisions. _Psychological Science, 9,_ 347-356. 
%
% * Ratcliff, R., & Rouder, J.N. (2000). A diffusion model account of masking
% in letter identification. _Journal of Experimental Psychology: Human
% Perception and Performance, 26,_ 127-140.
%
% * Ratcliff, R., Van Zandt, T., & McKoon, G. (1999). Connectionist and
% diffusion models of reaction time. _Psychological Review, 106,_ 261-300. 
%
% * Smith, P.L., & Vickers, D. (1988). The accumulator model to two-choice
% discrimination. _Journal of Mathematical Psychology, 32,_ 135-168. 
%
% * Strayer, D. L., & Kramer, A. F. (1994). Strategies and automaticity: I.
% Basic findings and conceptual framework. _Journal of Experimental
% Psychology: Learning, Memory, and Cognition, 20,_ 318-341.
%
% * Thapar, A., Ratcliff, R., & McKoon, G. (2003). A diffusion model analysis
% of the effects of aging on letter discrimination. _Psychology and Aging,
% 18,_ 415-429.
%
% * Townsend, J., & Ashby, F. (1983). _The stochastic modeling of elementary
% psychological processes._ Cambridge: Cambridge University Press. 
%
% * Tuerlinckx, F. (2004). The efficient computation of the distribution
% function of the diffusion process. _Behavior Research Methods,
% Instruments, & Computers, 36,_ 702-716.

%% Author of this file
% Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
% Part of the DMA Toolbox. Please read the End User License Agreement,
% contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
% See also http://ppw.kuleuven.be/okp/dmatoolbox.