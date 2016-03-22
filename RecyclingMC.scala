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
        println("mean squared displacement")
        println(r2.mean + "\t" + r2.error() + "\t" + r2.std())
    }
}

object RecyclingMC {
    def main(args: Array[String]) {
        var ho = new HarmonicOscillator(args(0).toDouble, args(1).toInt)
        ho.run()
  }
}
