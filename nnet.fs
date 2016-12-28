(* Author: Lance Zhou (Github: lazho)
 * Description: Neural network + genetic algorithm implementation
 *)

module nnet

type Neuron = float list

type NeuronLayer = Neuron list

type NeuronNetwork = NeuronLayer * NeuronLayer

let rand = System.Random()

(* Makes one neuron.
 * Params:
 *   ni Number of inputs
 *)
let makeNeuron (ni: int): Neuron =
  [for i in 1..(ni + 1) do yield rand.NextDouble()]

(* Makes one neuron layer.
 * Params:
 *   ni: Number of inputs
 *   nn: Number of neurons
 *)
let makeLayer (ni: int) (nn: int): NeuronLayer =
  let rec makeNeurons ni r =
    match r with
      | 0 -> []
      | _ -> (makeNeuron ni)::(makeNeurons ni (r-1))
  (makeNeurons ni nn)

(* Makes one network with one middle layer and one output layer.
 * Params:
 *   ni: Number of inputs
 *   no: Number of outputs
 *   nn: Number of middle neurons
 *)
let makeNetwork (ni: int) (no: int) (nn: int): NeuronNetwork =
  ((makeLayer ni nn), (makeLayer nn no))

(* Evaluates the output of a neuron given a set of inputs
 * Params:
 *   is: The set of inputs
 *   n:  The neuron
 *)
let evalNeuron is n =
  let rec evalHelper acc is n =
    match (is, n) with
      | ([], [wa]) -> (acc > wa)
      | (i::is, w::ws) ->
        if i then (evalHelper (acc+w) is ws)
        else (evalHelper acc is ws)
      | _ -> failwith "Incorrect input length."
  (evalHelper 0.0 is n)

(* Evaluates the outputs of a layer given a set of inputs
 * Params:
 *   is: The set of inputs
 *   l: The layer
 *)
let evalLayer is l =
  List.map (evalNeuron is) l

(* Evaluates the outputs of a two-layer network given a set
 * of inputs.
 * Params:
 *   is:  The set of inputs
 *   net: The network
 *)
let evalNetwork is net =
  let (ml, ol) = net in
    evalLayer (evalLayer is ml) ol

(* All the weights of a network flattened into one list, 
 * tupled with necessary info to reconstruct the chromosome
 * back into a network.
 *)
type Genome = float list * int * int * int

(* All the genomes coupled with their respective fitness
 * scores.
 *)
type GenePool = (Genome * int) list

(* Flatten a network into one list of weights, 
 * tupled with necessary info to reconstruct network
 * from chromosome.
 * Params:
 *   net: The given network
 *)
let netToGenome (net: NeuronNetwork): Genome =
  let (ml, ol) = net in
    match ml with
      | [] -> failwith "Empty middle layer"
      | n::ns ->
        (List.concat [(List.concat ml); (List.concat ol)],
          (n.Length - 1), (ol.Length), (ml.Length))

(* Lisp-like take implementation *)
let rec take l len =
  if (len = 0) then []
  else
    match l with
      | [] -> failwith "Ran out of list elements"
      | x::xs -> x::(take xs (len - 1))

(* Lisp-like drop implementation *)
let rec drop l len =
  if (len = 0) then l
  else
    match l with
      | [] -> failwith "Ran out of list elements"
      | x::xs -> (drop xs (len - 1))

(* Convenient take-drop helper *)
let splitList l len =
  ((take l len), (drop l len))

(* Reconstruct a network given its genome
 * Params:
 *   (ws, ni, no, nn): Deconstructed genome, see makeNetwork
 *)
let genomeToNet (ws, ni, no, nn): NeuronNetwork =
  let (ms, os) = (splitList ws ((ni + 1) * nn))
  let rec helper ns ni =
    match ns with
      | [] -> []
      | _ -> (take ns (ni + 1))::(helper (drop ns (ni + 1)) ni)
  ((helper ms ni), (helper os nn))

(* Mergesort part *)
let rec part l acc1 acc2 =
    match l with
        | [] -> (acc1, acc2)
        | x::xs -> (part xs (x::acc2) acc1)

(* Mergesort merge *)
let rec merge comp l r =
    match (l,r) with
        | ([],[]) -> []
        | (x::xs,[]) -> x::xs
        | ([],y::ys) -> y::ys
        | (x::xs,y::ys) ->
            if (comp x y) then x::(merge comp xs (y::ys))
            else y::(merge comp (x::xs) ys)

(* Mergesort main *)
let rec mergesort l comp =
    match l with
        | [] -> []
        | [x] -> [x]
        | x::xs ->
            let (l1,l2) = (part (x::xs) [] [])
            (merge comp (mergesort l1 comp) (mergesort l2 comp))

(* Use mergesort to sort gene pool for fittest genomes
 * Params:
 *   gp: The given gene pool
 *)
let sortByFittest (gp: GenePool): GenePool =
  (mergesort gp (fun (_, sx) (_, sy) -> sx > sy))

(* Construct a gene pool, given a list of genomes and the
 * fitness evaluation function.
 * Params:
 *   gs:  The list of genomes
 *   fit: The fitness function to evaluate the genomes with
 *)
let makeGenePool (gs: Genome list) (fit: NeuronNetwork -> int): GenePool =
  (sortByFittest (List.map (fun g -> (g, (fit (genomeToNet g)))) gs))

(* Helper function for crossover operation *)
let rec crossOverHelper adam eve rc =
  match (adam, eve) with
    | ([], []) -> []
    | (wa::was, we::wes) ->
      if ((rand.NextDouble()) > rc) then (we, wa)::(crossOverHelper was wes rc)
      else (wa, we)::(crossOverHelper was wes rc)
    | _ -> failwith "No cross species breeding allowed! (Adam, Eve not of same length)"

(* Helper function that takes ('a * 'a) list and makes 'a list list
 * by "unzipping" the tuples *)
let unzip l =
  let rec helperA l =
    match l with
      | [] -> []
      | (a, b)::xs -> a::(helperA xs)
  let rec helperB l =
    match l with
      | [] -> []
      | (a, b)::xs -> b::(helperB xs)
  [(helperA l); (helperB l)]

(* Crossover operation main.
 * Params:
 *   adam: First chromosome
 *   eve:  Second chromosome
 *   rc:   Crossover rate
 *)
let crossOver adam eve rc =
  (unzip (crossOverHelper adam eve rc))

(* Mutation operation
 * Params
 *   rm: Mutation rate
 *   g:  Chromosome to mutate
 *)
let rec mutate rm g =
  match g with
    | [] -> []
    | w::ws ->
      if ((rand.NextDouble()) > rm) then
        if ((rand.Next()) % 2 = 0) then (w + (rand.NextDouble()))::(mutate rm ws)
        else (w - (rand.NextDouble()))::(mutate rm ws)
      else w::(mutate rm ws)

(* Given a gene pool, take 2 fittest genomes and make 6 new based on them.
 * Return these 8 genomes in a new pool.
 * Params:
 *   fit: The function that evaluates fitness
 *   rc:  Crossover rate
 *   rm:  Mutation rate
 *   gp:  The current gene pool
 *)
let iterateGenePool fit rc rm (gp: GenePool) =
  let fittest = (take gp 2)
  match fittest with
    | (ag, _)::((eg, _)::[]) ->
      let ((adam, ai, ao, an), (eve, ei, eo, en)) = (ag, eg)
      (makeGenePool
        (List.map (fun ws -> (ws, ai, ao, an))
          (List.concat [[adam; eve]; 
            (List.map (mutate rm) 
              (List.concat [
                (crossOver adam eve rc);
                (crossOver adam eve rc);
                (crossOver adam eve rc)]))]))
        fit)
    | _ -> failwith "How did you get here?"

(* Helper function that turns a non-negative unsigned int into a list of bool
 * Used in the "recognise the number 4" test
 * Params:
 *   n:  The uint to turn into a bool list
 *   nb: The number of bools that are needed to represent n
 *)
let rec uintToBitList (n: int) (nb: int): bool list =
  if (nb = 0) then []
  else
    match (n / (pown 2 (nb - 1))) with
      | 0 -> false::(uintToBitList n (nb - 1))
      | 1 -> true::(uintToBitList (n - (pown 2 (nb - 1))) (nb - 1))
      | _ -> failwith "This isn't how this works at all" // There's no error checking.

(* Keep calling iterateGenePool until a minimum fitness score is reached by
 * the fittest genome in the pool.
 * Params:
 *   gc:  The current generation count. (Just set it to 0 when calling this)
 *   fit: The fitness eval function
 *   rc:  Crossover rate
 *   rm:  Mutation rate
 *   gp:  The current gene pool
 *   min: Minimum fitness score
 *)
let rec iterateUntilScoreReached gc fit rc rm gp min =
  match gp with
  | (g, s)::gs ->
    if (s < min) then
      (printfn "Generation #%i, highscore %i" gc s)
      (iterateUntilScoreReached (gc + 1) fit rc rm (iterateGenePool fit rc rm gp) min)
    else
      (genomeToNet g)
  | [] -> failwith "There is nothing in this gene pool."

(* Fitness eval function to check if a network can tell if a 3 by 5 array of pixels
 * is the number "4". The network has to take 15 inputs and have 1 output.
 * Params:
 *   net: The network being evaluated.
 *)
let testRecognise4 (net: NeuronNetwork): int =
  let mutable score = (pown 2 (3 * 5))
  for i in 0..((pown 2 (3 * 5)) - 1) do
    match (evalNetwork (uintToBitList i 15) net) with
      | [true] ->
        if (i <> 23497) then score <- score - 1
      | [false] ->
        if (i = 23497) then score <- score - 1
      | _ -> failwith "Too many net outputs"
  score
