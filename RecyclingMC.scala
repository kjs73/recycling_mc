import scala.math.exp
import scala.math.sqrt
import scala.util.Random

class MomentAcc() {
    var nr_samples: Int = 0
    var mean: Double = 0
    var mean2: Double = 0
    def add(inp: Double) {
        mean = (mean * nr_samples + inp) / (nr_samples + 1)
        mean2 = (mean2 * nr_samples + (inp * inp)) / (nr_samples + 1)
        nr_samples += 1
    }
    def error() : Double = {
        return sqrt((mean2 - mean * mean) / (nr_samples - 1))
    }
    def std() : Double = {
        return sqrt(mean2 - mean * mean);
    }
}

class HarmonicOscillator(kc: Double, nr_samplesc: Int) {
    val k: Double = kc
    val nr_samples: Int = nr_samplesc
    var x: Double = 0
    val beta: Double = 1
    var step_size: Double = 1
    val step_adapt_samples: Int = nr_samples / 2
    val equilibration_samples: Int = nr_samples / 2
    var r2 = new MomentAcc
    var acceptance = new MomentAcc
    var gen = Random
    gen.setSeed(42)
    def run() {
        choose_stepsize()
        take_steps(equilibration_samples)
        for (i <- 0 until nr_samples) {
            val accepted = take_step()
            record(accepted)
        }
        print_results()
    }
    def energy(state: Double) : Double = {
        return 0.5 * k * state * state
    }
    def take_step() : Boolean = {
        val x_new: Double = x + (0.5 - gen.nextDouble) * step_size
        if (gen.nextDouble < exp(-beta * (energy(x_new) - energy(x)))) {
            x = x_new
            return true
        }
        return false
    }
    def take_steps(N: Int) : Double = {
        var nr_accepted: Int = 0
        for (i <- 0 until N) {
            if (take_step()) {
                nr_accepted += 1
            }
        }
        return nr_accepted.toDouble / N
    }
    def record(accepted: Boolean) {
        if (accepted) {
            acceptance.add(1)
        }
        else {
            acceptance.add(0)
        }
        r2.add(x * x)
    }
    def choose_stepsize() {
        val adapt_batch: Int = 100
        val adapt_rounds: Int = step_adapt_samples / adapt_batch
        for (i <- 0 until adapt_rounds) {
            adapt_stepsize(adapt_batch, i)
        }
        take_steps(step_adapt_samples % adapt_batch) 
    }
    def adapt_stepsize(adapt_batch: Int, iteration: Int) {
        val acceptance_probability = take_steps(adapt_batch)
        val factor = (iteration + 1).toDouble / (iteration + 2)
        if (acceptance_probability < 0.2) {
            step_size *= factor
        }
        else if (acceptance_probability > 0.3) {
            step_size /= factor
        }
    }
    def print_results() {
        println("mean squared displacement, std")
        println(r2.mean + "\t" + r2.std())
        println("analytical result\n" + 1 / k)
        println("acceptance: " + acceptance.mean)
    }
}

class HarmonicOscillatorSymmetric(kc: Double, nr_samples: Int) extends HarmonicOscillator(kc, nr_samples) {
    override def take_step() : Boolean = {
        val x_new: Double = x + (0.5 - gen.nextDouble) * step_size
        if (gen.nextDouble < exp(-beta * energy(x_new)) / (exp(-beta * energy(x_new)) + exp(-beta * energy(x)))) {
            x = x_new
            return true
        }
        return false
    }
}

class HarmonicOscillatorRecycle(kc: Double, nr_samples: Int, nr_trial_movesc: Int) extends HarmonicOscillator(kc, nr_samples) {
    val nr_trial_moves: Int = nr_trial_movesc
    override def run() {
        choose_stepsize()
        take_steps(equilibration_samples)
        for (i <- 0 until nr_samples) {
            take_recycling_step()
        }
        print_results()
    }
    def take_recycling_step() {
        var trial_states: Array[Double] = new Array[Double](nr_trial_moves)
        var weights: Array[Double] = new Array[Double](nr_trial_moves)
        var r2_increment: Double = 0 
        var sum_weights: Double = 0
        for (i <- 0 until nr_trial_moves) {
            if (i > 0) {
                trial_states(i) = x + (0.5 - gen.nextDouble) * step_size
            }
            else {
                trial_states(0) = x
            }
            weights(i) = exp(-beta * energy(trial_states(i)))
            r2_increment += weights(i) * trial_states(i) * trial_states(i)
            sum_weights += weights(i)
        }
        r2_increment /= sum_weights
        r2.add(r2_increment)
        val threshold = gen.nextDouble * sum_weights
        var step_index: Int = 0
        var partial_sum = weights(0)
        while (partial_sum < threshold) {
            step_index += 1
            partial_sum += weights(step_index)
        }
        x = trial_states(step_index)
    }
    override def print_results() {
        println("mean squared displacement")
        println(r2.mean)
        println("analytical result\n" + 1 / k)
    }
}

object RecyclingMC {
    def main(args: Array[String]) {
        for (a <- args) {
            println(a)
        }
        if (args(0) == "metropolis") {
            new HarmonicOscillator(args(1).toDouble, args(2).toInt).run()
        }
        else if (args(0) == "symmetric") {
            new HarmonicOscillatorSymmetric(args(1).toDouble, args(2).toInt).run()
        }
        else if (args(0) == "recycle") {
            new HarmonicOscillatorRecycle(args(1).toDouble, args(2).toInt, args(3).toInt).run()
        }
        else {
            throw new IllegalArgumentException
        }
  }
}
