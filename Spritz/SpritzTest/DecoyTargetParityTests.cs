using NUnit.Framework;
using Proteomics;
using SpritzModifications;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace SpritzTest
{
    /// <summary>
    /// The target-decoy database has to hold one decoy per target. Issue #220 reported 90,584 decoys against
    /// 90,268 targets, meaning some decoys had no target counterpart and the decoy search space that false
    /// discovery rates are estimated against was larger than the target space it is meant to mirror.
    /// </summary>
    public class DecoyTargetParityTests
    {
        private static List<Protein> WithDecoys(params string[] sequences)
        {
            var targets = sequences.Select((s, i) => new Protein(s, "acc" + i)).ToList();

            return targets.Concat(DecoyProteinGenerator.GenerateDecoys(targets, DecoyType.Reverse)).ToList();
        }

        [Test]
        public void DropEntriesWithUnknownCTerminus_KeepsCompleteEntriesPaired()
        {
            var kept = ProteinAnnotation.DropEntriesWithUnknownCTerminus(WithDecoys("MPEPTIDEK", "MPEPTIDER", "SEQUENCEK"));

            Assert.That(kept.Count(p => !p.IsDecoy), Is.EqualTo(3));
            Assert.That(kept.Count(p => p.IsDecoy), Is.EqualTo(3));
        }

        /// <summary>
        /// The regression. An entry whose C-terminal residue is unknown is dropped, and its decoy has to go
        /// with it even though the decoy does not end in '?' - reversal moved that residue to the front.
        /// Covers both the initiator-methionine case and the case without one, since the two put the reversed
        /// C-terminus at different indices.
        /// </summary>
        [Test]
        public void DropEntriesWithUnknownCTerminus_DropsBothHalvesOfAnIncompleteEntry()
        {
            var kept = ProteinAnnotation.DropEntriesWithUnknownCTerminus(
                WithDecoys("MPEPTIDEK", "MPEPTIDE?", "SEQUENCE?"));

            Assert.That(kept.Count(p => !p.IsDecoy), Is.EqualTo(1), "only the complete target survives");
            Assert.That(kept.Count(p => p.IsDecoy), Is.EqualTo(1), "the incomplete entries must not leave decoys behind");
            Assert.That(kept.Select(p => p.BaseSequence), Has.None.Contains("?"),
                "no unknown residue should reach the database in either orientation");
        }

        /// <summary>
        /// Pins the assumption the filter rests on. If mzLib ever stops holding the initiator methionine in
        /// place, or reverses differently, this fails rather than silently unbalancing the database again.
        /// </summary>
        [Test]
        public void ReversalMovesTheUnknownResidueOffTheCTerminus()
        {
            var decoys = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { new Protein("MPEPTIDE?", "withM"), new Protein("SEQUENCE?", "withoutM") },
                DecoyType.Reverse);

            var withM = decoys.Single(p => p.Accession.Contains("withM"));
            var withoutM = decoys.Single(p => p.Accession.Contains("withoutM"));

            Assert.That(withM.BaseSequence.EndsWith('?'), Is.False);
            Assert.That(withM.BaseSequence[1], Is.EqualTo('?'), "held methionine puts the reversed C-terminus at index 1");
            Assert.That(withoutM.BaseSequence[0], Is.EqualTo('?'), "without a methionine it lands at index 0");
        }
    }
}
